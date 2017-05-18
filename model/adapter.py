import re

import aiohttp
import gnomic
import numpy as np
import json
from cobra import Metabolite, Reaction
from cobra.manipulation import find_gene_knockout_reactions
from cameo.data import metanetx

from model import logger
from model.settings import ID_MAPPER_API


def clean_bigg_id(string):
    return re.sub(r"bigg:|dsh", "", string)


async def query_identifiers(object_ids, db_from, db_to):
    """Call the id mapper service.

    :param object_ids: list of identifiers to query
    :param db_from: the source of the identifier, e.g. 'kegg'
    :param db_to: the destination type of the identifier, e.g. 'bigg'
    """
    query = json.dumps({'ids': object_ids, 'dbFrom': db_from, 'dbTo': db_to, 'type': 'Metabolite'})
    logger.info('query id mapper at {} with {}'.format(ID_MAPPER_API, str(query)))
    async with aiohttp.ClientSession() as session:
        async with session.post(ID_MAPPER_API, data=query) as r:
            assert r.status == 200
            result = await r.json()
            return result['ids']


def get_existing_metabolite(mnx_id, model, compartment):
    """Find compartment in the model by Metanetx id.

    Parameters
    ----------
    mnx_id : string
        Metanetx id
    model
        cobra model
    compartment : string
        f.e "_c"

    Returns
    -------
    Metabolite or None

    """
    if not mnx_id:
        return
    try:
        clean_id = clean_bigg_id(metanetx.mnx2bigg[mnx_id])
        return model.metabolites.get_by_id(clean_id + compartment)
    except KeyError:
        try:
            return model.metabolites.get_by_id(mnx_id + compartment)
        except KeyError:
            pass


def contains_carbon(metabolite):  # TODO: use method from Metabolite class when this change is merged
    if not metabolite.formula:
        raise ValueError("No formula for metabolite {}, it's unknown if there is carbon in it")
    return 'C' in metabolite.elements


def find_metabolite_info(met_id):
    """Find chemical formula of metabolite in metanetx.chem_prop dictionary

    Parameters
    ----------
    met_id : string
        string of format "<metabolite_id>_<compartment_id>", where <metabolite_id> is a BIGG id or a Metanetx id

    Returns
    -------
    pandas row or None

    """
    met_id = met_id[:-2]
    try:
        if met_id in metanetx.chem_prop.index:
            return metanetx.chem_prop.loc[met_id]
        return metanetx.chem_prop.loc[metanetx.all2mnx['bigg:' + met_id]]
    except KeyError:
        return None


def feature_id(feature):
    return feature.name if feature.name else feature.accession.identifier


class ModelModificationMixin(object):
    """
    Base model modification class, providing methods for adding adapter, exchange and demand reactions for new
    metabolites
    """
    model = None
    changes = None

    def create_exchange(self, metabolite):
        """For given metabolite A_<x> from compartment <x>, create:
        a) corresponding metabolite A_e from e compartment;
        b) adapter reaction A_<x> <--> A_e
        c) exchange reaction A_e -->

        Parameters
        ----------
        metabolite : Metabolite
            metabolite id in format <bigg_id>_<x>, f.e. Nacsertn_c

        Returns
        -------

        """
        extracellular_metabolite = Metabolite(metabolite.id[:2] + '_e', formula=metabolite.formula, compartment='e')
        extracellular_metabolite.added_by_model_adjustment = True
        self.add_adapter_reaction(metabolite, extracellular_metabolite)
        self.add_demand_reaction(extracellular_metabolite)

    def add_demand_reaction(self, metabolite):
        """Add demand reaction "A --> " for given metabolite A
         If reaction exists, log and pass

        Parameters
        ----------
        metabolite : basestring
            metabolite id in format <bigg_id>_<compartment_id>, f.e. Nacsertn_c

        Returns
        -------

        """
        try:
            logger.debug('Add demand reaction for metabolite: {}'.format(metabolite.id))
            demand_reaction = self.model.add_boundary(metabolite, type='demand')
            self.changes['added']['reactions'].add(demand_reaction)
            self.annotate_new_metabolites(demand_reaction)
        except ValueError:
            logger.debug('demand reaction exists for metabolite {}'.format(metabolite.id))

    def add_adapter_reaction(self, metabolite, existing_metabolite):
        """Add adapter reaction A <--> B for metabolites A and B

        Parameters
        ----------
        metabolite : Metabolite
            metabolite A
        existing_metabolite : Metabolite
            metabolite B

        Returns
        -------

        """
        try:
            adapter_reaction = Reaction(str('adapter_' + metabolite.id + '_' + existing_metabolite.id))
            adapter_reaction.lower_bound = -1000
            adapter_reaction.add_metabolites({metabolite: -1, existing_metabolite: 1})
            self.model.add_reactions([adapter_reaction])
            self.changes['added']['reactions'].add(adapter_reaction)
            self.annotate_new_metabolites(adapter_reaction)
            logger.debug('Adapter reaction added: {} <--> {}'.format(metabolite.id, existing_metabolite.id))
        except Exception:  # TODO: raise a reasonable exception on cobra side if the reaction exists
            logger.debug('Adapter reaction exists: {} <--> {}'.format(metabolite.id, existing_metabolite.id))

    def make_consumable(self, metabolite):
        """For metabolite in e compartment with existing exchange reaction, make it possible to consume metabolite
        by decreasing the lower bound of exchange reaction

        Parameters
        ----------
        metabolite : Metabolite
            metabolite from e compartment, f.e. melatn_e

        Returns
        -------

        """
        exchange_reaction = list(set(metabolite.reactions).intersection(self.model.exchanges))[0]
        if exchange_reaction.lower_bound >= 0:
            exchange_reaction.lower_bound = -1 if contains_carbon(metabolite) else -1000
        self.changes['added']['reactions'].add(exchange_reaction)
        self.annotate_new_metabolites(exchange_reaction)

    @staticmethod
    def annotate_new_metabolite(metabolite):
        """Find information about new metabolite in chem_prop dictionary and add it to model

        Parameters
        ----------
        metabolite : Metabolite
            new metabolite

        Returns
        -------

        """
        info = find_metabolite_info(metabolite.id)
        if info is not None:
            metabolite.formula = info['formula']
            metabolite.name = info['name']
            metabolite.annotation = info.to_dict()
        else:
            logger.debug('No formula for {}'.format(metabolite.id))

    def annotate_new_metabolites(self, reaction):
        """Annotate new metabolites with chem_prop information and keep track of them

        Parameters
        ----------
        reaction : Reaction
            reaction that is added to the model

        Returns
        -------

        """
        for metabolite in reaction.metabolites:
            if getattr(metabolite, 'added_by_model_adjustment', False):
                self.annotate_new_metabolite(metabolite)
                self.changes['added']['metabolites'].add(metabolite)

    def model_metabolite(self, metabolite_id, compartment='_e'):
        """Get metabolite associated with this model for a given entity

        Parameters
        ----------
        metabolite_id
            string of format <database>:<id>, f.e. chebi:12345
        compartment
            the compartment where to find the metabolite, e.g. _e for extracellular compartment
        Returns
        -------
        the model metabolite (or None if no matching found)
        """
        mnx_id = metanetx.all2mnx.get(metabolite_id)
        return get_existing_metabolite(mnx_id, self.model, compartment)


async def map_equation_to_model(equation, model_namespace, compartment=None):
    """Try to map given equation which contains KEGG ids to an equation with ids
    in the model's namespace.
    If metabolite does not exist in the BIGG database, use Metanetx id.
    If compartment is given, metabolites ids will have it as postfix.

    Example:
    Input: C00002 + C00033 <=> C00013 + C05993, compartment='_c'
    Output: atp_c + ac_c <=> ppi_c + MNXM4377_c

    :param equation: string
    :param compartment: f.e. "_c"
    :param model_namespace: The namespace the model uses for its metabolites
    :return:
    """
    array = equation.split()
    re_id = re.compile("^[A-Za-z][A-Za-z0-9]*$")
    to_map = [el for el in array if re.match(re_id, el)]
    logger.info('query ids: {}'.format(to_map))
    mapping = await query_identifiers(to_map, 'kegg', model_namespace)
    logger.info('response ids: {}'.format(mapping))
    result = []
    for i, el in enumerate(array):
        if not re.match(re_id, el):
            result.append(el)
        else:
            try:
                el = mapping[el][0]
            except KeyError:
                pass
            if compartment:
                el += compartment
            result.append(el)
    return ' '.join(result)


def full_genotype(genotype_changes):
    """Construct gnomic Genotype object from the list of strings with changes

    :param genotype_changes: list of changes, f.e. ['-tyrA::kanMX+', 'kanMX-']
    :return:
    """

    def chain(definitions, **kwargs):
        if not definitions:
            return gnomic.Genotype([])
        genotype = gnomic.Genotype.parse(definitions[0], **kwargs)
        for definition in definitions[1:]:
            genotype = gnomic.Genotype.parse(definition, parent=genotype, **kwargs)
        return genotype

    return chain(genotype_changes)


class GenotypeChangeModel(ModelModificationMixin):
    """
    Applies genotype change on cameo model
    """

    def __init__(self, model, genotype_changes, genes_to_reactions, namespace):
        """Initialize change model

        :param model: cameo model
        :param genotype_changes: gnomic.Genotype object
        :param genes_to_reactions: dictionary like {<gene name>: {<reaction id>: <reactions equation>, ...}, ...}
        """
        self.compartment = '_c'
        self.model = model
        self.namespace = namespace
        self.genes_to_reactions = genes_to_reactions
        self.changes = {
            'added': {'reactions': set(), 'metabolites': set()},  # reaction contain information about genes
            'removed': {'genes': set(), 'reactions': set()},
        }

    async def apply_changes(self, genotype_changes):
        """Apply genotype changes on initial model

        :param genotype_changes: gnomic.Genotype
        :return:
        """
        for change in genotype_changes.changes():
            if isinstance(change, gnomic.Mutation):
                await self.apply_mutation(change)
            if isinstance(change, gnomic.Plasmid):
                await self.add_plasmid(change)

    async def apply_mutation(self, mutation):
        """Apply mutations on initial model

        :param mutation: gnomic.Mutation
        :return:
        """
        if mutation.old:
            for feature in mutation.old.features():
                self.knockout_gene(feature)
        if mutation.new:
            for feature in mutation.new.features():
                await self.add_gene(feature)

    async def add_plasmid(self, plasmid):
        """Add plasmid features to the initial model.
        No plasmid instance in cameo, so changes are made in model genes and reactions directly

        :param plasmid: gnomic.Plasmid
        :return:
        """
        for feature in plasmid.features():
            await self.add_gene(feature)

    def knockout_gene(self, feature):
        """Perform gene knockout.
        Use feature name as gene name

        :param feature: gnomic.Feature
        :return:
        """
        genes = self.model.genes.query(feature.name, attribute="name")
        if genes:
            gene = genes[0]
            gene.knock_out()
            self.changes['removed']['genes'].add(gene)
            for reaction in find_gene_knockout_reactions(self.model, [gene]):
                self.changes['removed']['reactions'].add(reaction)
            logger.info('Gene knockout: {}'.format(gene.name))
        else:
            logger.info('Gene for knockout is not found: {}'.format(feature.name))

    async def add_gene(self, feature):
        """Perform gene insertion.
        Find all the reactions associated with this gene using KEGGClient and add them to the model

        :param feature: gnomic.Feature
        :return:
        """
        logger.info('Add gene: {}'.format(feature.name))
        identifier = feature_id(feature)
        if self.model.genes.query(identifier, attribute='name'):  # do not add if gene is already there
            logger.info('Gene {} exists in the model'.format(feature.name))
            return
        for reaction_id, equation in self.genes_to_reactions.get(identifier, {}).items():
            await self.add_reaction(reaction_id, equation, identifier)
        logger.info('Gene added: {}'.format(identifier))

    async def add_reaction(self, reaction_id, equation, gene_name):
        """Add new reaction by rn ID from equation, where metabolites defined by kegg ids.

        :param reaction_id: reaction rn ID
        :param equation: equation string, where metabolites are defined by kegg ids
        :param gene_name: gene name
        :return:
        """
        reaction = Reaction(reaction_id)
        self.model.add_reactions([reaction])
        equation = await map_equation_to_model(equation, self.namespace, self.compartment)
        logger.info('New reaction: {}'.format(equation))
        # TODO: adjust when build_reaction_from_string returns the new metabolites >>
        metabolites_before = {m.id for m in self.model.metabolites}
        reaction.build_reaction_from_string(equation)
        metabolites_after = {m.id for m in self.model.metabolites}
        # TODO: <<
        for added_metabolite in metabolites_after.difference(metabolites_before):
            new_metabolite = self.model.metabolites.get_by_id(added_metabolite)
            new_metabolite.added_by_model_adjustment = True
            self.create_exchange(new_metabolite)
        reaction.gene_reaction_rule = gene_name
        self.changes['added']['reactions'].add(reaction)


class MediumChangeModel(ModelModificationMixin):
    """
    Applies medium on cameo model
    """

    def __init__(self, model, medium):
        """
        Parameters
        ----------
        model
            cameo model
        medium
            list of dictionaries of format
            {'id': <compound id (<database>:<id>, f.e. chebi:12345)>, 'concentration': <compound concentration (float)>}
        """
        self.medium = medium
        self.model = model
        self.changes = {
            'added': {'reactions': set()},
        }
        self.apply_medium()

    def apply_medium(self):
        """For each metabolite in medium try to find corresponding metabolite in e compartment of the model.
        If metabolite is found, change the lower limit of the reaction to a negative number,
        so the model would be able to consume this compound.
        If metabolite is not found in e compartment, log and continue.
        """
        for compound in self.medium:
            existing_metabolite = self.model_metabolite(compound['id'], '_e')
            if not existing_metabolite:
                logger.info('No metabolite {}'.format(compound['id']))
                continue
            logger.info('Found metabolite {}'.format(compound['id']))
            self.make_consumable(existing_metabolite)


class MeasurementChangeModel(ModelModificationMixin):
    """
    Update constraints based on measured fluxes
    """

    def __init__(self, model, measurements):
        """

        Parameters
        ----------
        model
            cameo model
        measurements
            A list of dictionaries of format
            {'id': <metabolite id (<database>:<id>, f.e. chebi:12345)>, 'measurements': list(<measurement (float)>)}
        """
        self.measurements = measurements
        self.model = model
        self.changes = {
            'added': {'reactions': set()},
        }
        self.missing_in_model = []
        self.apply_flux_bounds()

    def apply_flux_bounds(self):
        """For each measured flux (production-rate / uptake-rate), constrain the model by setting upper and lower
        bound to either the max/min values of the measurements if less than three observations, otherwise to 97%
        normal distribution range i.e., mean +- 1.96 * stdev. """
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
            model_metabolite = self.model_metabolite(scalar['id'], '_e')
            if not model_metabolite:
                self.missing_in_model.append(scalar['id'])
                logger.info('Model is missing metabolite {}'.format(scalar['id']))
                return
            reaction = list(set(model_metabolite.reactions).intersection(self.model.exchanges))[0]
            reaction.bounds = lower_bound, upper_bound
            self.changes['added']['reactions'].add(reaction)
