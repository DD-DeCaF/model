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

import logging

import numpy as np
from cobra import Reaction
from cobra.io.dict import reaction_to_dict

from model.exceptions import NoIDMapping, PartNotFound
from model.ice_client import ICE
from model.modeling.cobra_helpers import get_unique_metabolite
from model.modeling.driven import minimize_distance
from model.modeling.gnomic_helpers import feature_id, full_genotype
from model.modeling.salts import MEDIUM_SALTS


logger = logging.getLogger(__name__)
ice = ICE()


def apply_medium(model, medium):
    """
    Apply a medium to a metabolic model.

    The medium is applied by adding salt mappings, resolving the corresponding
    exchange reaction for the resulting medium compounds, setting their uptake
    rate, and then applying that to the model.

    Parameters
    ----------
    model: cobra.Model
    medium: list(dict)
        The medium definition, a list of dicts matching the `MediumCompound`
        schema.

    Returns
    -------
    tuple (operations, errors)
        Operations is a list of model operations necessary to apply the medium
        to the model. See also the `Operations` schema.
        If errors is not an empty list, it was not possible to apply the medium.
        Errors then contains a list of string messages describing the
        problem(s).
    """
    operations = []
    errors = []

    # Detect salt compounds
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
    medium.extend([
        {'id': 'chebi:25517', 'name': 'nickel'},
        {'id': 'chebi:25368', 'name': 'molybdate'},
    ])

    # Create a map of exchange reactions and corresponding fluxes to apply to
    # the medium.
    medium_mapping = {}
    for compound in medium:
        try:
            metabolite = get_unique_metabolite(model, compound['id'], 'e', 'CHEBI')
        except NoIDMapping:
            errors.append(f"Cannot add medium compund '{compound['id']}' - metabolite not found in extracellular "
                          "compartment")
        else:
            exchange_reactions = metabolite.reactions.intersection(model.exchanges)
            if len(exchange_reactions) != 1:
                errors.append(f"Medium compound metabolite '{metabolite.id}' has {len(exchange_reactions)} exchange "
                              f"reactions in the model; expected 1")
                continue
            exchange_reaction = next(iter(exchange_reactions))

            # If someone already figured out the uptake rate for the compound, it's likely more accurate than our
            # assumptions, so keep it
            if exchange_reaction.id in model.medium:
                medium_mapping[exchange_reaction.id] = model.medium[exchange_reaction.id]
                continue

            if not metabolite.formula:
                logger.warning(f"No formula for metabolite '{metabolite.id}', cannot check if it is a carbon source")
                # If we don't know, it's most likely that the metabolite does not have a higher uptake rate than a
                # carbon source, so set the bound still to 10
                medium_mapping[exchange_reaction.id] = 10
            elif 'C' in metabolite.elements:
                # Limit the uptake rate for carbon sources to 10
                medium_mapping[exchange_reaction.id] = 10
            else:
                medium_mapping[exchange_reaction.id] = 1000

    # Apply the medium to the model, letting cobrapy deal with figuring out the correct bounds to change
    model.medium = medium_mapping

    # Add all exchange reactions to operations, to make sure any changed bounds is properly updated
    for reaction in model.exchanges:
        operations.append({
            'operation': 'modify',
            'type': 'reaction',
            'id': reaction.id,
            'data': reaction_to_dict(reaction),
        })

    return operations, errors


def apply_genotype(model, genotype_changes):
    """
    Apply genotype changes to a metabolic model.

    The genotype is first parsed by gnomic. Removed features (genes) are knocked
    out, while added features are added by looking up reaction definitions in
    ICE and adding those to the model.

    Parameters
    ----------
    model: cobra.Model
    genotype_changes: list(str)
        A list of genotype change strings parseable by gnomic. For example:
        ["-tyrA::kanMX+", "kanMX-"].

    Returns
    -------
    tuple (operations, errors)
        Operations is a list of model operations necessary to apply the medium
        to the model. See also the `Operations` schema.
        If errors is not an empty list, it was not possible to apply the medium.
        Errors then contains a list of string messages describing the
        problem(s).
    """
    operations = []
    errors = []

    genotype_changes = full_genotype(genotype_changes)

    # Apply feature operations
    for feature in genotype_changes.removed_features:
        feature_name = feature_id(feature)
        # Perform gene knockout. Use feature name as gene name
        try:
            gene = model.genes.query(lambda g: feature_name in (g.id, g.name))[0]
            gene.knock_out()
            operations.append({
                'operation': 'knockout',
                'type': 'gene',
                'id': gene.id,
            })
        except IndexError:
            logger.warning(f"Cannot knockout gene '{feature_name}', not found in the model")

    for feature in genotype_changes.added_features:
        feature_name = feature_id(feature)
        # Perform gene insertion.
        # Find all the reactions associated with this gene using KEGGClient and add them to the model
        if model.genes.query(lambda g: feature_name in (g.id, g.name)):
            # do not add if gene is already in the model
            logger.info(f"Not adding gene '{feature_name}', it already exists in the model")
            continue

        try:
            for reaction_id, equation in ice.get_reaction_equations(genotype=feature_name).items():
                logger.info(f"Adding reaction '{reaction_id}' for gene '{feature_name}'")

                # Add the reaction
                try:
                    reaction = Reaction(reaction_id)
                    reaction.gene_reaction_rule = feature_name
                    model.add_reaction(reaction)
                except ValueError:
                    # The reaction ID is already in the model
                    continue

                # Before building the reaction's metabolites, keep track of the existing ones to detect new
                # metabolites added to the model
                metabolites_before = set(model.metabolites)
                reaction.build_reaction_from_string(equation)
                new_metabolites = set(model.metabolites).difference(metabolites_before)

                operations.append({
                    'operation': 'add',
                    'type': 'reaction',
                    'data': reaction_to_dict(reaction),
                })

                # For all new metabolites, create a demand reaction so that it may leave the system
                for metabolite in new_metabolites:
                    demand_reaction = model.add_boundary(metabolite, type='demand')
                    operations.append({
                        'operation': 'add',
                        'type': 'reaction',
                        'data': reaction_to_dict(demand_reaction),
                    })
        except PartNotFound:
            logger.warning(f"Cannot add gene '{feature_name}', no gene-protein-reaction rules were found in ICE")

    return operations, errors


def apply_measurements(model, biomass_reaction, measurements):
    """
    Apply omics measurements to a metabolic model.

    For each measured flux (production-rate / uptake-rate), constrain the model
    by forcing their upper and lower bounds to the measured values.

    Parameters
    ----------
    model: cobra.Model
    biomass_reaction: str
        The id of the biomass reaction in the given model.
    measurements: list(dict)
        The measurements, a list of dicts matching the `Measurement` schema.

    Returns
    -------
    tuple (operations, errors)
        Operations is a list of model operations necessary to apply the
        measurements to the model. See also the `Operations` schema.
        If errors is not an empty list, it was not possible to apply the
        measurements.
        Errors then contains a list of string messages describing the
        problem(s).
    """
    operations = []
    errors = []

    # First, improve the fluxomics dataset by minimizing the distance to a feasible problem
    measurements = minimize_distance(model, biomass_reaction, measurements)

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
                metabolite = get_unique_metabolite(model, compound, 'e', 'CHEBI')
            except NoIDMapping:
                errors.append(f"Cannot find compound '{compound}' in the model")
            else:
                exchange_reactions = metabolite.reactions.intersection(model.exchanges)
                if len(exchange_reactions) != 1:
                    errors.append(f"Measured metabolite '{metabolite.id}' has {len(exchange_reactions)} exchange "
                                  f"reactions in the model; expected 1")
                    continue
                exchange_reaction = next(iter(exchange_reactions))

                # data is adjusted assuming a forward exchange reaction, x <-- (sign = -1), so if we instead actually
                # have <-- x, then multiply with -1
                direction = exchange_reaction.metabolites[metabolite]
                if direction > 0:
                    lower_bound, upper_bound = -1 * lower_bound, -1 * upper_bound
                exchange_reaction.bounds = lower_bound, upper_bound
                operations.append({
                    'operation': 'modify',
                    'type': 'reaction',
                    'id': exchange_reaction.id,
                    'data': reaction_to_dict(exchange_reaction),
                })
        elif scalar['type'] == 'reaction':
            try:
                reaction = model.reactions.get_by_id(scalar['id'])
            except KeyError:
                errors.append(f"Cannot find reaction '{scalar['id']}' in the model")
            else:
                reaction.bounds = lower_bound, upper_bound
                operations.append({
                    'operation': 'modify',
                    'type': 'reaction',
                    'id': reaction.id,
                    'data': reaction_to_dict(reaction),
                })
        elif scalar['type'] == 'growth-rate':
            reaction = model.reactions.get_by_id(biomass_reaction)
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
                errors.append(f"Cannot find reaction '{scalar['id']}' in the model")
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
    return operations, errors