import asyncio
import aiohttp
import math
from collections import Counter
import logging

from cameo import models, phenotypic_phase_plane
from cobra.io.dict import (reaction_to_dict, gene_to_dict, metabolite_to_dict,
                           reaction_from_dict, metabolite_from_dict)
import gnomic

from model.adapter import (GenotypeChangeModel, full_genotype, MediumChangeModel,
                           get_unique_metabolite, NoIDMapping, MeasurementChangeModel,
                           feature_id)
import model.constants as constants
from model.settings import ANNOTATIONS_API

LOGGER = logging.getLogger(__name__)

async def operate_on_reactions(model, reactions, key, apply_function, undo_function):
    if 'changes' not in model.notes:
        model.notes['changes'] = constants.get_empty_changes()
    current = model.notes['changes'][key]['reactions']
    applied = set([r['id'] for r in current])
    if key == 'removed':
        to_apply = set(reactions) - applied
        to_undo = [r for r in current
                   if r['id'] in (applied - set(reactions))]
    else:
        to_apply = [r for r in reactions if r['id'] not in applied]
        to_undo = [r for r in current if r['id'] in applied - set([i['id'] for i in reactions])]
    new_reactions = await apply_function(model, to_apply)
    removed = undo_function(model, to_undo)
    model.notes['changes'][key]['reactions'] = [r for r in current if r['id'] in applied - removed]
    model.notes['changes'][key]['reactions'].extend(new_reactions)
    return model


async def add_reaction_from_universal(model, reaction_id):
    reaction = models.metanetx_universal_model_bigg.reactions.get_by_id(reaction_id)
    reaction_string = reaction.build_reaction_string()
    adapter = GenotypeChangeModel(
        model,
        [],
        {None: {reaction_id: reaction_string}},
        model.notes['namespace'],
    )
    await adapter.map_metabolites()
    adapter.add_reaction(reaction_id, reaction_string, None)
    return collect_changes(adapter)


def collect_changes(modifications):
    model = modifications.model
    changes = model.notes.get('changes', constants.get_empty_changes())
    to_dict = dict(
        reactions=reaction_to_dict,
        genes=gene_to_dict,
        metabolites=metabolite_to_dict,
    )
    for action, by_action in modifications.changes.items():
        for entity, value in by_action.items():
            ids = set([i['id'] for i in changes[action][entity]])
            to_replace = set([i.id for i in value if i.id in ids])
            changes[action][entity] = [i for i in changes[action][entity] if i['id'] not in to_replace]
            changes[action][entity].extend([to_dict[entity](i) for i in value])
    model.notes['changes'] = changes
    return model


# TODO: eliminate when cobrapy is refactored
def build_string_from_metabolites(metabolites):
    """Generate a human readable reaction string"""

    def format(number):
        return "" if number == 1 else str(number).rstrip(".") + " "

    reactant_bits = []
    product_bits = []
    for name, coefficient in metabolites.items():
        if coefficient >= 0:
            product_bits.append(format(coefficient) + name)
        else:
            reactant_bits.append(format(abs(coefficient)) + name)

    return '{} <=> {}'.format(' + '.join(reactant_bits), ' + '.join(product_bits))


async def add_apply(model, to_apply):
    added = []
    before = {r['id'] for r in model.notes['changes']['added']['reactions']}
    for rn in to_apply:
        if rn['metabolites']:
            model = await add_reaction_from_string(
                model,
                rn['id'],
                build_string_from_metabolites(rn['metabolites'])
            )
        else:
            model = await add_reaction_from_universal(model, rn['id'])
    for reaction in model.notes['changes']['added']['reactions']:
        if reaction['id'] not in before:
            added.append(reaction_to_dict(model.reactions.get_by_id(reaction['id'])))
    return added


def count_metabolites(reactions):
    reactions = [r for r in reactions if not is_dummy(r['id'])]
    metabolites = [m for reaction in reactions for m in reaction['metabolites']]
    return dict(Counter(metabolites).most_common())


def add_undo(model, to_undo):
    all_metabolites_count = count_metabolites(model.notes['changes']['added']['reactions'])
    undo_metabolites_count = count_metabolites(to_undo)
    exchanges_to_keep = set()
    for k in all_metabolites_count:
        if all_metabolites_count[k] - undo_metabolites_count.get(k, 0) >= 1:
            if k.endswith('_c'):
                exchanges_to_keep.add('DM_' + k[:-2] + '_e')
                exchanges_to_keep.add('adapter_' + k + '_' + k[:-2] + '_e')
    final_undo = set([i['id'] for i in to_undo]) - exchanges_to_keep
    model.remove_reactions([model.reactions.get_by_id(i) for i in final_undo], remove_orphans=True)
    metabolites_after = {m.id for m in model.metabolites}
    model.notes['changes']['added']['metabolites'] = [m for m in model.notes['changes']['added']['metabolites'] if
                                                      m['id'] in metabolites_after]
    return final_undo


async def add_reaction_from_string(model, reaction_id, reaction_string):
    reaction_string = reaction_string.strip()
    adapter = GenotypeChangeModel(
        model,
        [],
        {None: {reaction_id: reaction_string}},
        model.notes['namespace']
    )
    await adapter.map_metabolites()
    adapter.add_reaction(reaction_id, reaction_string, None, '')
    return collect_changes(adapter)


async def apply_reactions_add(model, reactions_ids):
    return await operate_on_reactions(model, reactions_ids, 'added', add_apply, add_undo)


def knockout_undo(model, to_undo):
    for reaction in to_undo:
        if model.reactions.has_id(reaction['id']):
            model.reactions.get_by_id(reaction['id']).bounds = \
                reaction['lower_bound'], reaction['upper_bound']
    return {reaction['id'] for reaction in to_undo}


async def knockout_apply(model, to_apply):
    removed = []
    for r_id in to_apply:
        if model.reactions.has_id(r_id):
            removed.append(reaction_to_dict(model.reactions.get_by_id(r_id)))
            model.reactions.get_by_id(r_id).knock_out()
    return removed


async def apply_reactions_knockouts(model, reactions_ids):
    return await operate_on_reactions(model, reactions_ids, 'removed', knockout_apply, knockout_undo)


async def call_genes_to_reactions(genotype_features):
    """Make series of asynchronous calls for getting information about which reactions were added with new genes

    :param genotype_features: generator of new genes ids
    :return:
    """
    identifiers = list(new_features_identifiers(genotype_features))
    results = await asyncio.gather(*[
        query_genes_to_reaction(gene=identifier) for identifier in identifiers
    ])
    return {k: v for k, v in zip(identifiers, results)}


async def query_genes_to_reaction(gene):
    """Make asynchronous call for getting information about which reactions were added with the gene

    :param gene: gene identifier
    :return: reactions mapping {<rn ID>: <reaction string>}
    """
    LOGGER.info('Annotated gene at %s: %s', ANNOTATIONS_API, gene)
    async with aiohttp.ClientSession() as session:
        async with session.get(ANNOTATIONS_API, params={'geneId': gene}) as r:
            try:
                assert r.status == 200
            except Exception as ex:
                print('######## ANNOTATIONS_API params: ', ANNOTATIONS_API, {'geneId': gene})
                print('######## ANNOTATIONS_API exception:', await r.text())
                raise ex
            result = await r.json()
            return result.get('response', {})


async def apply_genotype_changes(model, genotype_changes):
    """Apply genotype changes to cameo model.

    :param model: cameo model
    :param genotype_changes: list of strings, f.e. ['-tyrA::kanMX+', 'kanMX-']
    :return:
    """
    LOGGER.info('Genotype changes %s', genotype_changes)
    genotype_features = full_genotype(genotype_changes)
    genes_to_reactions = await call_genes_to_reactions(genotype_features)
    LOGGER.info('Genes to reaction: %s', genes_to_reactions)
    change_model = GenotypeChangeModel(model, genotype_features, genes_to_reactions, model.notes['namespace'])
    await change_model.map_metabolites()
    change_model.apply_changes()
    return change_model


async def apply_medium_changes(model, medium):
    change_model = MediumChangeModel(model, medium)
    change_model.apply_medium()
    return change_model


def convert_mg_to_mmol(mg, formula_weight):
    return mg * (1 / formula_weight)


def convert_measurements_to_mmol(measurements, model):
    for value in measurements:
        if 'unit' not in value:
            continue
        if value['unit'] == 'mg':
            try:
                metabolite = get_unique_metabolite(model, value['id'], 'e', 'CHEBI')
            except NoIDMapping:
                continue
            else:
                value['measurements'] = [convert_mg_to_mmol(
                    point,
                    metabolite.formula_weight
                ) for point in value['measurements']]
                value['unit'] = 'mmol'
                LOGGER.info('Converted metabolite %s from mg to mmol', value['id'])
    return measurements


def fix_measurements_ids(measurements):
    IDS = {
        'chebi:42758': 'chebi:12965',
    }
    for value in measurements:
        if 'id' in value and value['id'] in IDS:
            # id is to be interpreted as compound_id
            value['id'] = IDS[value['id']]
    return measurements

async def apply_measurement_changes(model, measurements):
    measurement = convert_measurements_to_mmol(measurements, model)
    measurements = fix_measurements_ids(measurement)
    change_model = MeasurementChangeModel(model, measurements)
    change_model.apply_measurements()
    return change_model


async def modify_model(message, model):
    apply_functions = {
        constants.GENOTYPE_CHANGES: apply_genotype_changes,
        constants.MEDIUM: apply_medium_changes,
        constants.MEASUREMENTS: apply_measurement_changes,
    }

    for key in constants.REQUEST_KEYS:
        data = message.get(key, [])
        if data:
            modifications = await apply_functions[key](model, data)
            model = collect_changes(modifications)
    if constants.REACTIONS_ADD in message:
        model = await apply_reactions_add(model, message[constants.REACTIONS_ADD])
    if constants.REACTIONS_KNOCKOUT in message:
        model = await apply_reactions_knockouts(model, message[constants.REACTIONS_KNOCKOUT])
    return model


def restore_changes(model, changes):
    LOGGER.info('Changes to restore: %s', changes)
    model = apply_additions(model, changes['added'])
    model = apply_removals(model, changes['removed'])
    model = apply_measurements(model, changes['measured'])
    return model


def apply_additions(model, changes):
    model = add_metabolites(model, changes['metabolites'])
    model = add_reactions(model, changes['reactions'])
    return model


def apply_measurements(model, changes):
    for rxn in changes['reactions']:
        model.reactions.get_by_id(rxn['id']).bounds = rxn['lower_bound'], rxn['upper_bound']
    return model


def add_reactions(model, changes):
    reactions = [reaction_from_dict(r, model) for r in changes]
    to_add = []
    current = set()
    for reaction in reactions:
        if model.reactions.has_id(reaction.id):
            model.reactions.get_by_id(reaction.id).bounds = reaction.lower_bound, reaction.upper_bound
        elif reaction.id not in current:
            to_add.append(reaction)
            current.add(reaction.id)
    model.add_reactions(to_add)
    return model


def add_metabolites(model, changes):
    metabolites = [metabolite_from_dict(m) for m in changes]
    model.add_metabolites(metabolites)
    return model


def apply_removals(model, changes):
    model = remove_reactions(model, changes['reactions'])
    model = remove_genes(model, changes['genes'])
    return model


def remove_genes(model, changes):
    for gene in changes:
        gene = model.genes.query(gene['name'], attribute="name")[0]
        gene.knock_out()
    return model


def remove_reactions(model, changes):
    for r in changes:
        reaction = model.reactions.get_by_id(r['id'])
        reaction.knock_out()
    return model


def is_dummy(reaction_id):
    return reaction_id.startswith('DM') or reaction_id.startswith('adapter')


def product_reaction_variable(model, metabolite_id):
    """Find a medium exchange reaction in the model for the given metabolite id

    :param model: cameo model
    :param metabolite_id: string of format <database>:<id>, f.e. chebi:12345
    :return:
    """
    db_name, compound_id = metabolite_id.split(':')
    try:
        metabolite = get_unique_metabolite(model, compound_id, 'e', db_name)
        LOGGER.info('found model metabolite for %s %s', compound_id, db_name)
    except NoIDMapping:
        LOGGER.info('no model metabolite for %s %s', compound_id, db_name)
        return None
    exchange_reactions = list(set(metabolite.reactions).intersection(model.exchanges))
    if len(exchange_reactions) != 1:
        LOGGER.info('More than one exchange reaction in model %s for metabolite %s', model.id, metabolite_id)
    return exchange_reactions[0]


def phase_plane_to_dict(model, metabolite_id):
    """Return phenotypic phase plane results in format that is convenient for

    :param model: cameo model
    :param metabolite_id: string of format <database>:<id>, f.e. chebi:12345
    :return:
    """
    reaction = product_reaction_variable(model, metabolite_id)
    if not reaction:
        return {}
    ppp = phenotypic_phase_plane(model, [reaction]).data_frame.to_dict()
    result = {}
    for k, v in ppp.items():
        if k not in {'c_yield_lower_bound', 'c_yield_upper_bound',
                     'mass_yield_lower_bound', 'mass_yield_upper_bound'}:
            result[k] = [float(v[point]) if not math.isnan(float(v[point])) else None for point in sorted(v.keys())]
    return result


def new_features_identifiers(genotype_changes: gnomic.Genotype):
    """Extract identifiers for features which addition is defined in gnomic string

    :param genotype_changes: gnomic string with genotype changes
    :return:
    """
    for change in genotype_changes.changes():
        if isinstance(change, gnomic.Mutation):
            if change.new:
                for feature in change.new.features():
                    yield feature_id(feature)
        if isinstance(change, gnomic.Plasmid):
            for feature in change.features():
                yield feature_id(feature)
