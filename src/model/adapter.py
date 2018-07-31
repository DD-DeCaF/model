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
import json
import logging
import re
from collections import defaultdict

import gnomic
import networkx as nx
import numpy as np
import pandas as pd
from cameo.data import metanetx
from cobra import Metabolite, Reaction
from cobra.manipulation import find_gene_knockout_reactions

from model import storage
from model.driven import adjust_fluxes2model
from model.exceptions import NoIDMapping
from model.gnomic_helpers import feature_id, insert_feature
from model.id_mapper_client import query_identifiers
from model.model_helpers import get_unique_metabolite, strip_compartment
from model.salts import MEDIUM_SALTS


logger = logging.getLogger(__name__)


def adapt_from_medium(model_meta, medium):
    return []


def adapt_from_genotype(model_meta, genotype):
    return []


def adapt_from_measurements(model_meta, measurements):
    return []


class ModelModificationMixin(object):
    """
    Base model modification class, providing methods for adding adapter, exchange and demand reactions for new
    metabolites
    """
    model = None
    changes = None
    metabolite_mapping = None

    def create_exchange(self, metabolite):
        """For given metabolite A_<x> from compartment <x>, create:
        a) corresponding metabolite A_e from e compartment;
        b) adapter reaction A_<x> <--> A_e
        c) exchange reaction A_e -->

        :param metabolite: cobra.Metabolite, metabolite id in format <bigg_id>_<x>, f.e. Nacsertn_c
        """
        extracellular_metabolite = Metabolite(re.sub('_[cpm]$', '_e', metabolite.id),
                                              formula=metabolite.formula, compartment='e')
        extracellular_metabolite.added_by_model_adjustment = True
        self.add_adapter_reaction(metabolite, extracellular_metabolite)
        self.add_demand_reaction(extracellular_metabolite)

    def add_demand_reaction(self, metabolite):
        """Add demand reaction "A --> " for given metabolite A
         If reaction exists, log and pass

        :param metabolite: basestring, metabolite id in format <bigg_id>_<compartment_id>, f.e. Nacsertn_c
        """
        try:
            logger.debug('Add demand reaction for metabolite: %s', metabolite.id)
            demand_reaction = self.model.add_boundary(metabolite, type='demand')
            self.changes['added']['reactions'].add(demand_reaction)
            self.annotate_new_metabolites(demand_reaction)
        except ValueError:
            logger.debug('demand reaction exists for metabolite %s', metabolite.id)

    def add_adapter_reaction(self, metabolite, existing_metabolite):
        """Add adapter reaction A <--> B for metabolites A and B

        :param metabolite: cobra.Metabolite, metabolite A
        :param existing_metabolite: cobra.Metabolite, metabolite B
        """
        try:
            adapter_reaction = Reaction(str('adapter_' + metabolite.id + '_' + existing_metabolite.id))
            adapter_reaction.lower_bound = -1000
            adapter_reaction.add_metabolites({metabolite: -1, existing_metabolite: 1})
            self.model.add_reactions([adapter_reaction])
            self.changes['added']['reactions'].add(adapter_reaction)
            self.annotate_new_metabolites(adapter_reaction)
            logger.debug('Adapter reaction added: %s <--> %s', metabolite.id, existing_metabolite.id)
        except Exception:  # TODO: raise a reasonable exception on cobra side if the reaction exists
            logger.debug('Adapter reaction exists: %s <--> %s', metabolite.id, existing_metabolite.id)

    def make_consumable(self, metabolite, key='measured'):
        """For metabolite in e compartment with existing exchange reaction, make it possible to consume metabolite
        by decreasing the lower bound of exchange reaction

        :param metabolite: cobra.Metabolite, metabolite from e compartment, f.e. melatn_e
        :param key: where to save the reactions, default is 'measured'
        """
        def contains_carbon(metabolite):  # TODO: use method from Metabolite class when this change is merged
            if not metabolite.formula:
                logger.warning(f"No formula for metabolite {metabolite.id}, it's unknown if there is carbon in it. "
                               "Assuming that there is no carbon")
                return False
            return 'C' in metabolite.elements

        exchange_reaction = list(set(metabolite.reactions).intersection(self.model.exchanges))[0]
        if exchange_reaction.lower_bound >= 0:
            exchange_reaction.lower_bound = -10 if contains_carbon(metabolite) else -1000
        self.changes[key]['reactions'].add(exchange_reaction)
        self.annotate_new_metabolites(exchange_reaction)

    def annotate_new_metabolite(self, metabolite):
        """Find information about new metabolite in chem_prop dictionary and add it to model

        :param metabolite: cobra.Metabolite, new metabolite
        """

        def find_key_for_id(met_id, metabolite_mapping):
            """Find mapped key for a given metabolite id.

             metabolite_mapping could e.g. {'MNXM89795': {'mnx': ['MNXM89795'], 'bigg': ['udpgal'], 'chebi': [
             '18307', '13487', '13495'], ...} and the metabolite.id 'udpgal'. We want to use this mapping to lookup
             the chebi from udpgal, this function will fetch the corresponding primary key 'MNXM89795'.
            """
            keys = [key for key, key_map in metabolite_mapping.items() if met_id in
                    set(itertools.chain(*key_map.values()))]
            if len(keys) == 0:
                return None
            elif len(keys) == 1:
                return keys[0]
            else:
                raise KeyError('multiply mapping keys for {}'.format(met_id))

        def find_metabolite_info(met_id):
            """Find chemical formula of metabolite in metanetx.chem_prop dictionary

            :param met_id: string, string of format "<metabolite_id>_<compartment_id>", where <metabolite_id> is a BIGG id or
            a Metanetx id
            :returns: pandas row or None
            """
            met_id = strip_compartment(met_id)
            try:
                if met_id in metanetx.chem_prop.index:
                    return metanetx.chem_prop.loc[met_id]
                return metanetx.chem_prop.loc[metanetx.all2mnx['bigg:' + met_id]]
            except KeyError:
                return None

        info = find_metabolite_info(metabolite.id)
        if info is not None:
            metabolite.formula = info['formula']
            metabolite.name = info['name']
            metabolite.annotation = info.to_dict()
        else:
            logger.debug('no formula for %s', metabolite.id)
        if self.metabolite_mapping is not None:
            mapped_id = find_key_for_id(strip_compartment(metabolite.id), self.metabolite_mapping)
            if mapped_id is not None:
                metabolite.annotation['CHEBI'] = [
                    f'{prefix}:{i}' if not i.startswith(f"{prefix}:") else i
                    for i in self.metabolite_mapping[mapped_id].get('chebi', [])
                ]
                metabolite.annotation['metanetx.chemical'] = self.metabolite_mapping[mapped_id].get('mnx', [])
                metabolite.annotation['bigg.metabolite'] = self.metabolite_mapping[mapped_id].get('bigg', [])
            else:
                logger.debug('no cross-references for %s', metabolite.id)
        logger.info('new annotation: %s', metabolite.annotation)

    def annotate_new_metabolites(self, reaction):
        """Annotate new metabolites with chem_prop information and keep track of them

        :param reaction: cobra.Reaction, reaction that is added to the model
        """
        for metabolite in reaction.metabolites:
            if getattr(metabolite, 'added_by_model_adjustment', False):
                self.annotate_new_metabolite(metabolite)
                self.changes['added']['metabolites'].add(metabolite)


class GenotypeChangeModel(ModelModificationMixin):
    """
    Applies genotype change on cameo model
    """

    def __init__(self, model, genotype_changes, genes_to_reactions, namespace):
        """Initialize change model

        :param model: cameo model
        :param genotype_changes: gnomic.Genotype object
        :param genes_to_reactions: dictionary like {<gene name>: {<reaction id>: <reactions equation>, ...}, ...}
        :param namespace: the namespace for the model's identifiers
        """
        self.compartment = '_c'
        self.model = model
        self.genes_to_reactions = genes_to_reactions
        self.genotype_changes = genotype_changes
        self.namespace = namespace
        self.metabolite_mapping = defaultdict(dict)
        self.metabolite_identifiers = defaultdict(list)

        # TODO: the metabolite identifiers in genes_to_reactions should be prefixed by the db_name
        re_id = re.compile("^[A-Za-z0-9_-]+$")
        re_id_map = {'mnx': re.compile(r'MNXM[\d]+'),
                     'kegg': re.compile(r'\bC[0-9]{5}\b'),
                     'bigg': re.compile(r'[a-z0-9_]+', re.IGNORECASE)}
        for reactions in self.genes_to_reactions.values():
            for equation in reactions.values():
                for element in equation.split():
                    element = strip_compartment(element)
                    if not element.isdigit() and re.match(re_id, element):
                        for db_key, pattern in re_id_map.items():
                            if re.match(pattern, element):
                                self.metabolite_identifiers[db_key].append(element)
                                break
        self.changes = {
            'added': {'reactions': set(), 'metabolites': set()},  # reaction contain information about genes
            'removed': {'genes': set(), 'reactions': set()},
        }

    def map_metabolites(self):
        logger.info("Metabolites to map: %s", self.metabolite_identifiers)
        from_db_key = list(self.metabolite_identifiers.keys())
        to_db_key = list({self.namespace, 'mnx', 'chebi'})
        queries = [(db_from, db_to) for db_from, db_to in itertools.product(from_db_key, to_db_key) if db_from != db_to]
        mappings = [query_identifiers(self.metabolite_identifiers[db_from], db_from, db_to)
                    for db_from, db_to in queries]
        for query, mapping in zip(queries, mappings):
            db_from, db_to = query
            for met_id in self.metabolite_identifiers[db_from]:
                self.metabolite_mapping[met_id][db_from] = [met_id]
                if met_id in mapping:
                    self.metabolite_mapping[met_id][db_to] = mapping[met_id]
        logger.info('Using metabolite mapping: %s', self.metabolite_mapping)

    def apply_changes(self):
        """Apply genotype changes on initial model

        :return:
        """
        def detect_mutations(remove, add):
            """Check if deletion-addition combination is in fact mutation and should be skipped.
            If the feature has DELETE_GENE suffix, consider this gene deleted
            """
            DELETE_GENE = 'delta8bp'  # the gene mutation which disables its functions
            drop_remove, drop_add = [], []
            for old_feature in remove:
                for new_feature in add:
                    if new_feature.startswith(old_feature):
                        logger.info('Gene %s is replaced with %s', old_feature, new_feature)
                        drop_add.append(new_feature)
                        if DELETE_GENE not in new_feature:
                            drop_remove.append(old_feature)
            for name in drop_remove:
                remove.pop(name)
            for name in drop_add:
                add.pop(name)

        to_remove = {}
        to_add = {}

        for change in self.genotype_changes.changes():
            if isinstance(change, gnomic.Mutation):
                old = change.old.features() if change.old else []
                for feature in old:
                    insert_feature(feature, to_remove, to_add)
                new = change.new.features() if change.new else []
                for feature in new:
                    insert_feature(feature, to_add, to_remove)
            if isinstance(change, gnomic.Plasmid):
                for feature in change.features():
                    insert_feature(feature, to_add, to_remove)

        detect_mutations(to_remove, to_add)

        for k, v in to_remove.items():
            self.knockout_gene(v)
        for k, v in to_add.items():
            self.add_gene(v)

    def knockout_gene(self, feature):
        """Perform gene knockout.
        Use feature name as gene name

        :param feature: gnomic.Feature
        :return:
        """
        genes = self.model.genes.query(feature.name, attribute="name")
        if feature.name in self.model.genes:
            genes += self.model.genes.get_by_any(feature.name)
        if genes:
            gene = genes[0]
            gene.knock_out()
            self.changes['removed']['genes'].add(gene)
            for reaction in find_gene_knockout_reactions(self.model, [gene]):
                self.changes['removed']['reactions'].add(reaction)
            logger.info('Gene knockout: %s', gene.name)
        else:
            logger.info('Gene for knockout is not found: %s', feature.name)

    def add_gene(self, feature):
        """Perform gene insertion.
        Find all the reactions associated with this gene using KEGGClient and add them to the model

        :param feature: gnomic.Feature
        :return:
        """
        logger.info('Add gene: %s', feature.name)
        identifier = feature_id(feature)
        if self.model.genes.query(identifier, attribute='name'):  # do not add if gene is already there
            logger.info('Gene %s exists in the model', feature.name)
            return
        for reaction_id, equation in self.genes_to_reactions.get(identifier, {}).items():
            self.add_reaction(reaction_id, equation, identifier)
        logger.info('Gene added: %s', identifier)

    def add_reaction(self, reaction_id, equation, gene_name, compartment=None):
        """Add new reaction by rn ID from equation, where metabolites defined by kegg ids.

        :param reaction_id: reaction rn ID
        :param equation: equation string, where metabolites are defined by kegg ids
        :param gene_name: gene name
        :return:
        """
        def map_equation_to_model(equation, metabolite_mapping, native_namespace='bigg', compartment=''):
            """Map given equation with some metabolite identifiers to those used in the model's namespace.

            If compartment is given, metabolites ids will have it as postfix.

            Example:
            Input: C00002 + C00033 <=> C00013 + C05993, compartment='_c', {'C00002': 'atp_c', 'C00033': 'ac_c', ...}
            Output: atp_c + ac_c <=> ppi_c + MNXM4377_c

            :param equation: string
            :param compartment: f.e. "_c"
            :param metabolite_mapping: A dict that maps metabolite identifiers to those used in the model
            :param native_namespace: string, the model's namespace to prefer for the identifiers
            :return:
            """

            def map_metabolite(met_id):
                for namespace in [native_namespace, 'mnx']:
                    if namespace in metabolite_mapping[met_id]:
                        return metabolite_mapping[met_id][namespace][0]
                return met_id

            array = equation.split()
            re_id = re.compile("^[A-Za-z0-9_-]+$")
            result = []
            for el in array:
                if el.isdigit() or not re.match(re_id, el):
                    result.append(el)
                else:
                    result.append(map_metabolite(el) + compartment)
            return ' '.join(result)

        if compartment is None:
            compartment = self.compartment
        if self.model.reactions.has_id(reaction_id):
            return
        reaction = Reaction(reaction_id)
        self.model.add_reactions([reaction])
        logger.info('New reaction to add: %s', equation)
        equation = map_equation_to_model(equation, self.metabolite_mapping, self.model.notes['namespace'], compartment)
        logger.info('New adjusted reaction: %s', equation)
        # TODO: adjust when build_reaction_from_string returns the new metabolites >>
        metabolites_before = {m.id for m in self.model.metabolites}
        reaction.build_reaction_from_string(equation)
        metabolites_after = {m.id for m in self.model.metabolites}
        # TODO: <<
        for added_metabolite in metabolites_after.difference(metabolites_before):
            new_metabolite = self.model.metabolites.get_by_id(added_metabolite)
            new_metabolite.added_by_model_adjustment = True
            self.create_exchange(new_metabolite)
        if gene_name:
            reaction.gene_reaction_rule = gene_name
        self.changes['added']['reactions'].add(reaction)
        self.annotate_new_metabolites(reaction)


class MediumChangeModel(ModelModificationMixin):
    """
    Applies medium on cameo model
    """
    TRACE_METALS = [
        {'id': 'chebi:25517', 'name': 'nickel'},
        {'id': 'chebi:25368', 'name': 'molybdate'},
    ]

    def __init__(self, model, medium):
        """
        :param model: cobra.Model
        :param medium: list of dictionaries of format
            {'id': <compound id (<database>:<id>, f.e. chebi:12345)>, 'concentration': <compound concentration (float)>}
        """
        self.medium = medium
        self.model = model
        self.changes = {
            'measured-medium': {'reactions': set()},
        }

    def detect_salt_compounds(self):
        result = []
        chebi_ids = [c['id'] for c in self.medium]
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
                            result.append({
                                'id': 'chebi:' + compound,
                            })
        return result

    def apply_medium(self):
        """For each metabolite in medium try to find corresponding metabolite in e compartment of the model.
        If metabolite is found, change the lower limit of the reaction to a negative number,
        so the model would be able to consume this compound.
        If metabolite is not found in e compartment, log and continue.
        """
        old_medium = []
        for r in self.model.medium:
            reaction = self.model.reactions.get_by_id(r)
            reaction.lower_bound = 0
            old_medium.append(reaction)
        self.medium.extend(self.detect_salt_compounds())
        self.medium.extend(self.TRACE_METALS)
        for compound in self.medium:
            try:
                existing_metabolite = get_unique_metabolite(self.model, compound['id'], 'e', 'CHEBI')
            except NoIDMapping:
                logger.info('No metabolite %s in external compartment', compound['id'])
                continue
            logger.info('Found metabolite %s', compound['id'])
            self.make_consumable(existing_metabolite, key='measured-medium')
        for reaction in old_medium:
            if reaction.id not in self.changes['measured-medium']['reactions']:
                self.changes['measured-medium']['reactions'].add(reaction)


class MeasurementChangeModel(ModelModificationMixin):
    """
    Update constraints based on measured fluxes
    """

    def __init__(self, model, measurements):
        """
        :param model: cobra.Model
        :param measurements:
            A list of dictionaries of format
            {'id': <metabolite id (<database>:<id>, f.e. chebi:12345)>, 'measurements': list(<measurement (float)>)}
        """
        self.measurements = measurements
        self.model = model
        self.changes = {
            'measured': {'reactions': set()},
            'added': {'reactions': set()},
            'measured-missing': {'reactions': set()}
        }
        self.missing_in_model = []

    def has_transport(self, metabolite_id, direction):
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
        for c in self.model.compartments:
            m_id = metabolite_id + '_' + c
            if self.model.metabolites.has_id(m_id):
                metabolites.append(self.model.metabolites.get_by_id(m_id))
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

    def allow_transport(self, metabolite_e, direction):
        """
        If transport between cytosol and extracellular space in desired direction
        does not already exist, create a helper transport reaction.
        """
        met_id = metabolite_e.id.replace('_e', '')
        metabolite_c = self.model.metabolites.get_by_id(met_id + '_c')
        if self.has_transport(met_id, direction):
            return
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
        self.model.add_reactions([transport_reaction])
        self.changes['added']['reactions'].add(transport_reaction)

    def reaction_for_compound(self, compound, lower_bound, upper_bound):
        try:
            model_metabolite = get_unique_metabolite(self.model, compound, 'e', 'CHEBI')
        except NoIDMapping:
            self.missing_in_model.append(compound)
            logger.info('Model is missing metabolite %s, cannot apply measurement', compound)
            return
        possible_reactions = list(set(model_metabolite.reactions).intersection(self.model.exchanges))
        if len(possible_reactions) > 1:
            logger.warn('using first of %s', ', '.join([r.id for r in possible_reactions]))
        reaction = possible_reactions[0]
        # data is adjusted assuming a forward exchange reaction, x <-- (sign = -1), so if we instead actually
        # have <-- x, then multiply with -1
        direction = reaction.metabolites[model_metabolite]
        self.allow_transport(model_metabolite, lower_bound)
        if direction > 0:
            lower_bound, upper_bound = -1 * lower_bound, -1 * upper_bound
        reaction.bounds = lower_bound, upper_bound
        return reaction

    def reaction_for_scalar(self, scalar, lower_bound, upper_bound):
        reaction = None
        if scalar['db_name'] != 'bigg.reaction':
            raise NotImplementedError('only supporting bigg reaction identifiers not %s' % scalar['db_name'])
        try:
            reaction = self.model.reactions.get_by_id(scalar['id'])
        except KeyError:
            self.changes['measured-missing']['reactions'].add(Reaction(scalar['id']))
        else:
            reaction.bounds = lower_bound, upper_bound
        return reaction

    def protein_exchange_reaction_for_scalar(self, scalar, upper_bound):
        reaction = None
        if scalar['db_name'] != 'uniprot':
            raise NotImplementedError('only supporting uniprot protein identifiers not %s' % scalar['db_name'])
        try:
            def query_fun(rxn):
                xrefs = rxn.annotation.get(scalar['db_name'], [])
                xrefs = xrefs if isinstance(xrefs, list) else [xrefs]
                return scalar['id'] in xrefs

            reaction = self.model.reactions.query(query_fun)[0]
        except (IndexError, KeyError):
            self.changes['measured-missing']['reactions'].add(Reaction('prot_{}_exchange'.format(scalar['id'])))
        else:
            reaction.bounds = 0, upper_bound
        return reaction

    def minimize_distance(self):
        """Replaces fluxomics measurements with the minimized distance"""
        index = []
        observations = []
        uncertainties = []

        for measure in [m for m in self.measurements if m['type'] == 'reaction']:
            index.append(measure['id'])
            observations.append(np.nanmean(measure['measurements']))
            if len(measure['measurements']) >= 3:
                uncertainties.append(np.nanstd(measure['measurements'], ddof=1))
            else:
                uncertainties.append(1)

        try:
            # Include growth rate measurement
            growth_rate = next(m for m in self.measurements if m['type'] == 'growth-rate')
            index.append(storage.get(self.model.id).growth_rate_reaction)
            observations.append(np.nanmean(growth_rate['measurements']))
            uncertainties.append(np.nanstd(growth_rate['measurements'], ddof=1)
                                 if len(growth_rate['measurements']) >= 3 else 1)
        except StopIteration:
            pass

        observations = pd.Series(index=index, data=observations)
        uncertainties = pd.Series(index=index, data=uncertainties)

        solution = adjust_fluxes2model(self.model, observations, uncertainties)
        for reaction, minimized_distance in solution.fluxes.iteritems():
            for measurement in self.measurements:
                if (measurement['type'] == 'growth-rate' and reaction == storage.get(self.model.id).growth_rate_reaction
                        or reaction == measurement.get('id')):
                    measurement['measurements'] = [minimized_distance]

    def apply_measurements(self):
        """For each measured flux (production-rate / uptake-rate), constrain the model by setting upper and lower
        bound to either the max/min values of the measurements if less than three observations, otherwise to 97%
        normal distribution range i.e., mean +- 1.96 * stdev. """
        def next_measured_reaction(exchange_reaction):
            """Find the reaction flux of which would be fully defined by the corresponding exchange reaction.

            :param exchange_reaction: cobra.Reaction
            :return: cobra.Reaction or None
            """
            met = list(exchange_reaction.metabolites.keys())[0]
            connected_reactions = [r for r in met.reactions if r != exchange_reaction]
            if len(connected_reactions) == 0:
                logger.warning('No reactions except exchange is connected to %s', met)
                return
            if len(connected_reactions) > 1:
                return
            next_reaction = connected_reactions[0]
            metabolites = [m for m in next_reaction.metabolites if m != met]
            if len(metabolites) > 1:
                return
            met2 = metabolites[0]
            if next_reaction.metabolites[met] * next_reaction.metabolites[met2] > 0:
                return
            return next_reaction

        for scalar in self.measurements:
            scalar_data = np.array([v for v in scalar['measurements'] if not np.isnan(v)])
            if len(scalar_data) > 2:
                upper_bound = float(np.mean(scalar_data) + 1.96 * np.std(scalar_data, ddof=1))
                lower_bound = float(np.mean(scalar_data) - 1.96 * np.std(scalar_data, ddof=1))
            elif len(scalar_data) > 0:
                upper_bound = float(np.max(scalar_data))
                lower_bound = float(np.min(scalar_data))
            else:
                continue
            reaction = None
            if scalar['type'] == 'compound':
                reaction = self.reaction_for_compound(scalar['id'], lower_bound, upper_bound)
            elif scalar['type'] == 'reaction':
                reaction = self.reaction_for_scalar(scalar, lower_bound, upper_bound)
            elif scalar['type'] == 'growth-rate':
                reaction = self.model.reactions.get_by_id(storage.get(self.model.id).growth_rate_reaction)
                reaction.bounds = lower_bound, upper_bound
            elif scalar['type'] == 'protein' and scalar['mode'] == 'quantitative':
                reaction = self.protein_exchange_reaction_for_scalar(scalar, upper_bound)
            else:
                logger.info('scalar for measured type %s not supported', scalar['type'])
            if reaction:
                self.changes['measured']['reactions'].add(reaction)
                if scalar['type'] == 'compound':
                    next_measured = next_measured_reaction(reaction)
                    if next_measured:
                        self.changes['measured']['reactions'].add(next_measured)
