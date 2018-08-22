# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import itertools
import logging

import gnomic
import networkx as nx
import numpy as np
from cobra import Metabolite, Reaction
from cobra.io.dict import reaction_to_dict

from model import storage
from model.driven import minimize_distance
from model.exceptions import NoIDMapping, PartNotFound
from model.gnomic_helpers import full_genotype, insert_feature
from model.ice_client import ICE
from model.model_helpers import get_unique_metabolite
from model.salts import MEDIUM_SALTS


logger = logging.getLogger(__name__)
ice = ICE()


def adapt_from_medium(model, medium):
    """
    Returns a list of operations to apply the given medium to the given model

    :param model: cobra.Model
    :param medium: list of dictionaries of format
        {'id': <compound id (<database>:<id>, f.e. chebi:12345)>, 'concentration': <compound concentration (float)>}
    """

    operations = []

    for reaction_id in model.medium:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.lower_bound = 0
        operations.append({
            'operation': 'modify',
            'type': 'reaction',
            'id': reaction.id,
            'data': reaction_to_dict(reaction),
        })

    # Detect salt compounds
    # NOTES(Ali): could be simplified by generating a different datastructre from `update_salts
    chebi_ids = [c['id'] for c in medium]
    for compound in chebi_ids:
        chebi_id = compound.replace('chebi:', '')
        if chebi_id in MEDIUM_SALTS:
            compounds = MEDIUM_SALTS[chebi_id]
            n_not_found = len([i for i in compounds if not i])
            logger.info('Metabolite %s can be splitted up to %s', chebi_id, str(compounds))
            if n_not_found:
                logger.info('For %s %d mappings of %d is not found', chebi_id, n_not_found, len(compounds))
            for array in compounds:
                for compound in array:
                    if compound:
                        medium.append({'id': 'chebi:' + compound})

    # Add trace metals
    # NOTES(Ali): why these hardcoded metals?
    medium.extend([
        {'id': 'chebi:25517', 'name': 'nickel'},
        {'id': 'chebi:25368', 'name': 'molybdate'},
    ])

    # Make all metabolites in the medium consumable by setting the exchange reactions lower bound to a negative number
    for compound in medium:
        try:
            metabolite = get_unique_metabolite(model, compound['id'], 'e', 'CHEBI')
        except NoIDMapping:
            # NOTES(Ali): why not add the metabolite? maybe we don't know which reactions to bind it to? do we need to?
            # NOTES(Ali): if it's an ionic compound then both it and its salts will be attempted added to the medium
            logger.info(f"Cannot add medium compund '{compound['id']}' - metabolite not found in in extracellular "
                        "compartment")
            continue
        else:
            exchange_reaction = list(set(metabolite.reactions).intersection(model.exchanges))[0]
            # NOTES(Ali): Why only if >= 0?
            if exchange_reaction.lower_bound >= 0:
                if not metabolite.formula:
                    logger.warning(f"No formula for metabolite {metabolite.id}, it's unknown if there is carbon in it. "
                                   "Assuming that there is no carbon")
                if not metabolite.formula or 'C' not in metabolite.elements:
                    exchange_reaction.lower_bound = -1000
                else:
                    exchange_reaction.lower_bound = -10
            operations.append({
                'operation': 'modify',
                'type': 'reaction',
                'id': exchange_reaction.id,
                'data': reaction_to_dict(exchange_reaction),
            })
            # NOTES(Ali): annotation has been removed here, hope that's okay?
    return operations


def adapt_from_genotype(model, genotype_changes):
    """
    Return a list of operations to apply to a model based on the given genotype changes.

    :param model: cobra.Model
    :param genotype_changes: list of genotype change strings, f.e. ['-tyrA::kanMX+', 'kanMX-']
    """
    # NOTES(Ali): metabolite mappings are completely removed now, what are the consequences?
    # NOTES(Ali): review gene duplication logic below (`insert_feature` vs duplication check)

    operations = []

    features_to_remove = {}
    features_to_add = {}

    for change in full_genotype(genotype_changes).changes():
        if isinstance(change, gnomic.Mutation):
            old = change.old.features() if change.old else []
            for feature in old:
                insert_feature(feature, features_to_remove, features_to_add)
            new = change.new.features() if change.new else []
            for feature in new:
                insert_feature(feature, features_to_add, features_to_remove)
        if isinstance(change, gnomic.Plasmid):
            for feature in change.features():
                insert_feature(feature, features_to_add, features_to_remove)

    # Check if deletion-addition combination is in fact mutation and should be skipped.
    for old_feature in features_to_remove:
        for new_feature in features_to_add:
            if new_feature.startswith(old_feature):
                logger.info(f"Dropping new feature '{new_feature}' because it was already removed")
                del features_to_add[new_feature]

    # Perform gene knockout. Use feature name as gene name
    for feature in features_to_remove.values():
        try:
            gene = model.genes.query(lambda g: feature.name in (g.id, g.name))[0]
            gene.knock_out()
            operations.append({
                'operation': 'remove',
                'type': 'gene',
                'id': gene.id,
            })
            # NOTES(Ali): removed `find_gene_knockout_reactions` added to changes - adding the gene removal operation
            # NOTES(Ali): directly instead
        except IndexError:
            logger.warning(f"Cannot knockout gene '{feature.name}', not found in the model")

    # Perform gene insertion.
    # Find all the reactions associated with this gene using KEGGClient and add them to the model
    for feature in features_to_add.values():
        feature_id = feature.name or feature.accession.identifier
        if model.genes.query(lambda g: feature.name in (g.id, g.name)):
            # do not add if gene is already in the model
            logger.info(f"Not adding gene '{feature.name}', it already exists in the model")
            continue

        try:
            # :param equation: equation string, where metabolites are defined by kegg ids
            for reaction_id, equation in ice.get_reaction_equations(genotype=feature_id).items():
                logger.info(f"Adding reaction '{reaction_id}' for gene '{feature_id}'")

                # Add the reaction
                try:
                    reaction = Reaction(reaction_id)
                    reaction.gene_reaction_rule = feature_id
                    model.add_reaction(reaction)
                except ValueError:
                    # The reaction ID is already in the model
                    continue

                # Before building the reaction's metabolites, keep track of the existing ones to detect new metabolites
                # added to the model
                metabolites_before = set(model.metabolites)
                reaction.build_reaction_from_string(equation)
                new_metabolites = set(model.metabolites).difference(metabolites_before)

                operations.append({
                    'operation': 'add',
                    'type': 'reaction',
                    'id': reaction.id,
                    'data': reaction_to_dict(reaction),
                })

                # NOTES(Ali): not mapping equation ids to model namespace now

                # For all new metabolites, create a demand reaction so that it may leave the system
                for metabolite in new_metabolites:
                    demand_reaction = model.add_boundary(metabolite, type='demand')
                    operations.append({
                        'operation': 'add',
                        'type': 'reaction',
                        'id': demand_reaction.id,
                        'data': reaction_to_dict(demand_reaction),
                    })
        except PartNotFound:
            logger.warning(f"Cannot add gene '{feature_id}', no gene-protein-reaction rules were found in ICE")

    return operations


def adapt_from_measurements(model, measurements):
    """
    For each measured flux (production-rate / uptake-rate), constrain the model by forcing their upper and lower bounds
    to the measured values.

    :param model: cobra.Model
    :param measurements:
        A list of dictionaries of format
        {'id': <metabolite id (<database>:<id>, f.e. chebi:12345)>, 'measurements': list(<measurement (float)>)}
    """
    operations = []

    # First, improve the fluxomics dataset by minimizing the distance to a feasible problem
    measurements = minimize_distance(model, measurements)

    for scalar in measurements:
        # If there are three or more observations, use the 97% normal distribution range, i.e., mean +- 1.96.
        # Otherwise, just use the max/min values of the measurements.
        scalar_data = np.array([v for v in scalar['measurements'] if not np.isnan(v)])
        if len(scalar_data) > 2:
            upper_bound = float(np.mean(scalar_data) + 1.96 * np.std(scalar_data, ddof=1))
            lower_bound = float(np.mean(scalar_data) - 1.96 * np.std(scalar_data, ddof=1))
        elif len(scalar_data) > 0:
            upper_bound = float(np.max(scalar_data))
            lower_bound = float(np.min(scalar_data))
        else:
            continue

        if scalar['type'] == 'compound':
            try:
                compound = scalar['id']
                model_metabolite = get_unique_metabolite(model, compound, 'e', 'CHEBI')
            except NoIDMapping:
                # NOTES(Ali): how to deal with measurement not found in the model?
                continue
            else:
                possible_reactions = list(set(model_metabolite.reactions).intersection(model.exchanges))
                if len(possible_reactions) > 1:
                    logger.warn('using first of %s', ', '.join([r.id for r in possible_reactions]))
                reaction = possible_reactions[0]

                operations.extend(_allow_transport(model, model_metabolite, lower_bound))

                # data is adjusted assuming a forward exchange reaction, x <-- (sign = -1), so if we instead actually
                # have <-- x, then multiply with -1
                direction = reaction.metabolites[model_metabolite]
                if direction > 0:
                    lower_bound, upper_bound = -1 * lower_bound, -1 * upper_bound
                reaction.bounds = lower_bound, upper_bound
                operations.append({
                    'operation': 'modify',
                    'type': 'reaction',
                    'id': reaction.id,
                    'data': reaction_to_dict(reaction),
                })
        elif scalar['type'] == 'reaction':
            try:
                reaction = model.reactions.get_by_id(scalar['id'])
            except KeyError:
                # NOTES(Ali): how to deal with measurement not found in the model?
                continue
            else:
                reaction.bounds = lower_bound, upper_bound
                operations.append({
                    'operation': 'modify',
                    'type': 'reaction',
                    'id': reaction.id,
                    'data': reaction_to_dict(reaction),
                })
        elif scalar['type'] == 'growth-rate':
            reaction = model.reactions.get_by_id(storage.get(model.id).growth_rate_reaction)
            reaction.bounds = lower_bound, upper_bound
            operations.append({
                'operation': 'modify',
                'type': 'reaction',
                'id': reaction.id,
                'data': reaction_to_dict(reaction),
            })
        elif scalar['type'] == 'protein' and scalar['mode'] == 'quantitative':
            try:
                def query_fun(rxn):
                    xrefs = rxn.annotation.get(scalar['db_name'], [])
                    xrefs = xrefs if isinstance(xrefs, list) else [xrefs]
                    return scalar['id'] in xrefs

                reaction = model.reactions.query(query_fun)[0]
            except (IndexError, KeyError):
                # NOTES(Ali): how to deal with measurement not found in the model?
                continue
            else:
                reaction.bounds = 0, upper_bound
                operations.append({
                    'operation': 'modify',
                    'type': 'reaction',
                    'id': reaction.id,
                    'data': reaction_to_dict(reaction),
                })
        else:
            raise NotImplementedError(f"Measurement type '{scalar['type']}' is not supported")

    return operations


def _has_transport(model, metabolite_id, direction):
    """
    Check if transport between cytosol and extracellular space in desired direction exists.
    Use the graph of the transport reaction between the same metabolite
    in different compartments.
    """
    G = nx.DiGraph()
    m_c, m_e = metabolite_id + '_c', metabolite_id + '_e'
    G.add_node(m_c)
    G.add_node(m_e)
    metabolites = []
    for c in model.compartments:
        m_id = metabolite_id + '_' + c
        if model.metabolites.has_id(m_id):
            metabolites.append(model.metabolites.get_by_id(m_id))
    for a, b in itertools.combinations(metabolites, 2):
        for reaction in set(a.reactions).intersection(set(b.reactions)):
            a_coef, b_coef = reaction.metabolites[a], reaction.metabolites[b]
            if a_coef * b_coef < 0:
                l_bound, r_bound = reaction.bounds
                if l_bound < 0 and r_bound > 0:
                    G.add_edge(a.id, b.id)
                    G.add_edge(b.id, a.id)
                else:
                    if b_coef > 0:
                        pair = a.id, b.id
                    else:
                        pair = b.id, a.id
                    if l_bound < 0:
                        G.add_edge(pair[1], pair[0])
                    else:
                        G.add_edge(*pair)
    if direction > 0:
        return nx.has_path(G, m_c, m_e)
    return nx.has_path(G, m_e, m_c)


def _allow_transport(model, metabolite_e, direction):
    """
    If transport between cytosol and extracellular space in desired direction
    does not already exist, create a helper transport reaction.
    """
    # NOTES(Ali): can we do the same as for the `create_exchange` logic here?
    # NOTES(Ali): different phrasing: for a measured compound, what is really required here?
    met_id = metabolite_e.id.replace('_e', '')
    metabolite_c = model.metabolites.get_by_id(met_id + '_c')
    if _has_transport(model, met_id, direction):
        return []
    if direction > 0:
        m_from, m_to = metabolite_c, metabolite_e
    else:
        m_from, m_to = metabolite_e, metabolite_c
    logger.info('Transport reaction for %s is not found, '
                'creating a transport reaction from %s to %s',
                met_id, m_from, m_to)
    transport_reaction = Reaction('adapter_' + met_id)
    transport_reaction.bounds = -1000, 1000
    transport_reaction.add_metabolites(
        {m_from: -1, m_to: 1}
    )
    model.add_reaction(transport_reaction)
    return [{
        'operation': 'add',
        'type': 'reaction',
        'id': transport_reaction.id,
        'data': reaction_to_dict(transport_reaction),
    }]
