import asyncio
import aiohttp_cors
import aiohttp
import gnomic
import json
import hashlib
import time
import aioredis
import os
import re
from copy import deepcopy
from collections import namedtuple
from itertools import chain
from functools import lru_cache
from aiohttp import web, WSMsgType
from cameo.data import metanetx
from cameo import phenotypic_phase_plane
from cameo import models, load_model as cameo_load_model
from cameo.flux_analysis import room, lmoma, moma
from cobra.flux_analysis import pfba, flux_variability_analysis
from cobra.exceptions import OptimizationError
from cameo.util import ProblemCache
from cobra.io import read_sbml_model
from cobra.io.dict import (model_to_dict, reaction_to_dict, reaction_from_dict, gene_to_dict,
                           metabolite_to_dict, metabolite_from_dict)
from model.adapter import (get_existing_metabolite, GenotypeChangeModel, MediumChangeModel,
                           MeasurementChangeModel, full_genotype, feature_id, query_identifiers)
from model import logger
from model.settings import ANNOTATIONS_API

SPECIES_TO_MODEL = {
    'ECOLX': ['iJO1366', 'e_coli_core'],
    'YEAST': ['iMM904', 'ecYeast7'],
    'CRIGR': ['iMM1415'],
    'CORGT': ['iNJ661'],
    'PSEPU': ['iJN746'],
}

MODELS = frozenset(chain.from_iterable(models for _, models in SPECIES_TO_MODEL.items()))
LOCAL_MODELS = frozenset(['ecYeast7'])

MODEL_NAMESPACE = {
    'iJO1366': 'bigg',
    'iMM904': 'bigg',
    'iMM1415': 'bigg',
    'iNJ661': 'bigg',
    'iJN746': 'bigg',
    'e_coli_core': 'bigg',
    'ecYeast7': 'yeast7'
}

MODEL_GROWTH_RATE = {
    'iJO1366': 'BIOMASS_Ec_iJO1366_core_53p95M',
    'iMM904': 'BIOMASS_SC5_notrace',
    'iMM1415': 'BIOMASS_mm_1_no_glygln',
    'iNJ661': 'BIOMASS_Mtb_9_60atp',
    'iJN746': 'BIOMASS_KT_TEMP',
    'e_coli_core': 'BIOMASS_Ecoli_core_w_GAM',
    'ecYeast7': 'r_2111'
}


def pfba_fva(model, reactions=None):
    return flux_variability_analysis(
        model,
        fraction_of_optimum=1,
        pfba_factor=1.05,
        reactions_list=reactions
    )


METHODS = {
    'fba': lambda model: model.optimize(),
    'pfba': pfba,
    'fva': flux_variability_analysis,
    'pfba-fva': pfba_fva,
    'room': room,
    'moma': moma,
    'lmoma': lmoma,
}

GENOTYPE_CHANGES = 'genotype-changes'
MEDIUM = 'medium'
MEASUREMENTS = 'measurements'
SIMULATION_METHOD = 'simulation-method'
MAP = 'map'
REACTIONS_KNOCKOUT = 'reactions-knockout'
REACTIONS_ADD = 'reactions-add'
MODEL = 'model'
FLUXES = 'fluxes'
GROWTH_RATE = 'growth-rate'
TMY = 'tmy'
OBJECTIVES = 'objectives'
REQUEST_ID = 'request-id'
REMOVED_REACTIONS = 'removed-reactions'
ADDED_REACTIONS = 'added-reactions'
MISSING_MEASURED_REACTIONS = 'missing-measured-reactions'

EMPTY_CHANGES = {
    'added': {
        'reactions': [],
        'metabolites': [],
    },
    'removed': {
        'genes': [],
        'reactions': [],
    },
    'measured': {
        'genes': [],
        'reactions': []
    },
    'measured-missing': {
        'genes': [],
        'reactions': []
    }
}

MAPS_DIR = 'maps'


def generate_map_dictionary():
    """Generate model-maps lookup depending on the folder structure

    :return: dict
    """
    result = {}
    for path, _, files in os.walk(MAPS_DIR):
        if files:
            result[path.replace(MAPS_DIR + '/', '')] = \
                sorted([re.match(r".*\.(.+)\..*", f).group(1) for f in files])
    return result


MAP_DICTIONARY = generate_map_dictionary()


async def redis_client():
    return await aioredis.create_redis((os.environ['REDIS_ADDR'], os.environ['REDIS_PORT']),
                                       loop=asyncio.get_event_loop())


def load_model(model_id):
    if model_id in LOCAL_MODELS:
        sbml_file = os.path.join(os.path.dirname(__file__), 'data', model_id + '.sbml.gz')
        model = read_sbml_model(sbml_file)
    else:
        model = cameo_load_model(model_id)
    model.notes['namespace'] = MODEL_NAMESPACE[model_id]
    return model


class Models(object):
    MODELS = {
        v: load_model(v) for v in MODELS
        }
    print('Models are ready')


async def find_changes_in_db(model_id):
    dumped_changes = None
    redis = await redis_client()
    if await redis.exists(model_id):
        dumped_changes = (await redis.get(model_id)).decode('utf-8')
    redis.close()
    return dumped_changes


async def restore_from_db(model_id):
    t = time.time()
    changes = await find_changes_in_db(model_id)
    if not changes:
        return None
    model = model_from_changes(changes)
    t = time.time() - t
    logger.info('Model with db key {} is ready in {} sec'.format(model_id, t))
    return model


@lru_cache(maxsize=2 ** 6)
def model_from_changes(changes):
    changes = json.loads(changes)
    model = find_in_memory(changes['model']).copy()
    model.notes = deepcopy(model.notes)
    model = restore_changes(model, changes['changes'])
    model.notes['changes'] = changes['changes']
    return model


def find_in_memory(model_id):
    return Models.MODELS.get(model_id)


async def restore_model(model_id):
    """Try to restore model by model id.
    NOTE: if model is found in memory, the original model is returned - to modify, make a copy

    :param model_id: str
    :return: Cameo model or None
    """
    model = find_in_memory(model_id)
    if model:
        logger.info('Wild type model with id {} is found'.format(model_id))
        return model
    model = await restore_from_db(model_id)
    if model:
        logger.info('Model with id {} found in database'.format(model_id))
        return model
    logger.info('No model with id {}'.format(model_id))
    return None


def key_from_model_info(wild_type_id, message):
    """Generate hash string from model information which will later be used as db key

    :param wild_type_id: str
    :param message: dict
    :return: str
    """
    d = {k: message.get(k, []) for k in {GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS}}
    d['model_id'] = wild_type_id
    return hashlib.sha224(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()


async def save_changes_to_db(model, wild_type_id, message):
    """Store model in cache database

    :param model: Cameo model
    :param wild_type_id: str
    :param message: dict
    :return: str (cache database key)
    """
    db_key = key_from_model_info(wild_type_id, message)
    redis = await redis_client()
    await redis.set(db_key,
                    json.dumps({'model': model.id, 'changes': model.notes.get('changes', deepcopy(EMPTY_CHANGES))}))
    redis.close()
    logger.info('Model created on the base of {} with message {} saved as {}'.format(wild_type_id, message, db_key))
    return db_key


class NoIDMapping(Exception):
    def __init__(self, metabolite_id):
        self.value = metabolite_id

    def __str__(self):
        return 'No Metanetx mapping for metabolite {}'.format(self.value)


async def existing_metabolite(model, metabolite_id):
    """Find metabolite in _e compartment of the given model.

    :param model: cameo model
    :param metabolite_id: string of format <database>:<id>, f.e. chebi:12345
    :return:
    """
    try:
        mnx_id = metanetx.all2mnx[metabolite_id]
    except KeyError:
        raise NoIDMapping(metabolite_id)
    elements_map = {
        'MNXM161568': 'MNXM322',  # thiamine
        'MNXM7043': 'MNXM41',  # glucose
    }
    if mnx_id in elements_map:
        mnx_id = elements_map[mnx_id]
    return await get_existing_metabolite(mnx_id, model, '_e')


async def product_reaction_variable(model, metabolite_id):
    """Find a medium exchange reaction in the model for the given metabolite id

    :param model: cameo model
    :param metabolite_id: string of format <database>:<id>, f.e. chebi:12345
    :return:
    """
    metabolite = await existing_metabolite(model, metabolite_id)
    if not metabolite:
        return None
    exchange_reactions = list(set(metabolite.reactions).intersection(model.exchanges))
    if len(exchange_reactions) != 1:
        logger.info('More than one exchange reaction in model {} for metabolite {}'.format(model.id, metabolite_id))
    return exchange_reactions[0]


async def phase_plane_to_dict(model, metabolite_id):
    """Return phenotypic phase plane results in format that is convenient for

    :param model: cameo model
    :param metabolite_id: string of format <database>:<id>, f.e. chebi:12345
    :return:
    """
    model.solver = 'glpk'
    reaction = await product_reaction_variable(model, metabolite_id)
    if not reaction:
        return {}
    ppp = phenotypic_phase_plane(model, [reaction]).data_frame.to_dict()
    result = {}
    for k, v in ppp.items():
        if k not in {'c_yield_lower_bound', 'c_yield_upper_bound',
                     'mass_yield_lower_bound', 'mass_yield_upper_bound'}:
            result[k] = [float(v[point]) for point in sorted(v.keys())]
    model.solver = 'cplex'
    return result


async def apply_genotype_changes(model, genotype_changes):
    """Apply genotype changes to cameo model.

    :param model: cameo model
    :param genotype_changes: list of strings, f.e. ['-tyrA::kanMX+', 'kanMX-']
    :return:
    """
    logger.info('Genotype changes {}'.format(genotype_changes))
    genotype_features = full_genotype(genotype_changes)
    genes_to_reactions = await call_genes_to_reactions(genotype_features)
    logger.info('Genes to reaction: {}'.format(genes_to_reactions))
    change_model = GenotypeChangeModel(model, genotype_features, genes_to_reactions, model.notes['namespace'])
    await change_model.map_metabolites()
    change_model.apply_changes()
    return change_model


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
    logger.info('Annotated gene at {}: {}'.format(ANNOTATIONS_API, gene))
    async with aiohttp.ClientSession() as session:
        async with session.get(ANNOTATIONS_API, params={'geneId': gene}) as r:
            assert r.status == 200
            result = await r.json()
            return result.get('response', {})


def convert_mg_to_mmol(mg, formula_weight):
    return mg * (1 / formula_weight)


async def apply_measurement_changes(model, measurements):
    measurement = await convert_measurements_to_mmol(measurements, model)
    measurements = fix_measurements_ids(measurement)
    change_model = MeasurementChangeModel(model, measurements)
    await change_model.apply_flux_bounds()
    return change_model


async def convert_measurements_to_mmol(measurements, model):
    for value in measurements:
        if 'unit' not in value:
            continue
        if value['unit'] == 'mg':
            metabolite = await existing_metabolite(model, value['id'])
            if metabolite:
                value['measurements'] = [convert_mg_to_mmol(
                    point,
                    metabolite.formula_weight
                ) for point in value['measurements']]
                value['unit'] = 'mmol'
                logger.info('Converted metabolite {} from mg to mmol'.format(value['id']))
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


async def apply_medium_changes(model, medium):
    change_model = MediumChangeModel(model, medium)
    await change_model.apply_medium()
    return change_model


ReactionKnockouts = namedtuple('ReactionKnockouts', ['model', 'changes'])


async def operate_on_reactions(model, reactions_ids, key, apply_function, undo_function):
    if 'changes' not in model.notes:
        model.notes['changes'] = deepcopy(EMPTY_CHANGES)
    current = model.notes['changes'][key]['reactions']
    applied = set([r['id'] for r in current])
    to_apply = set(reactions_ids) - applied
    to_undo = [r for r in current
               if r['id'] in (applied - set(reactions_ids))]
    new_reactions = await apply_function(model, to_apply)
    removed = undo_function(model, to_undo)
    model.notes['changes'][key]['reactions'] = [r for r in current if r['id'] in applied - removed]
    model.notes['changes'][key]['reactions'].extend(new_reactions)
    return model


async def apply_reactions_knockouts(model, reactions_ids):
    return await operate_on_reactions(model, reactions_ids, 'removed', knockout_apply, knockout_undo)


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


async def add_reaction_from_universal(model, reaction_id):
    reaction = models.metanetx_universal_model_bigg.reactions.get_by_id(reaction_id)
    reaction_string = reaction.build_reaction_string()
    adapter = GenotypeChangeModel(
        model,
        [],
        {None: {reaction_id: reaction_string}},
        model.notes['namespace']
    )
    await adapter.map_metabolites(from_namespace='mnx', template=r'MNXM[\d]+')
    adapter.add_reaction(reaction_id, reaction_string, None)
    return collect_changes(adapter)


async def apply_reactions_add(model, reactions_ids):
    return await operate_on_reactions(model, reactions_ids, 'added', add_apply, add_undo)


def is_dummy(reaction_id):
    return reaction_id.startswith('DM') or reaction_id.startswith('adapter')


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
    model.notes['changes']['added']['metabolites'] = [m for m in model.notes['changes']['added']['metabolites'] if m['id'] in metabolites_after]
    return final_undo


def count_metabolites(reactions):
    metabolites_count = {}
    for reaction in reactions:
        if not is_dummy(reaction['id']):
            for metabolite in reaction['metabolites']:
                metabolites_count[metabolite] = metabolites_count.get(metabolite, 0) + 1
    return metabolites_count


async def add_apply(model, to_apply):
    added = []
    for r_id in to_apply:
        model = await add_reaction_from_universal(model, r_id)
        for reaction in model.notes['changes']['added']['reactions']:
            added.append(reaction_to_dict(model.reactions.get_by_id(reaction['id'])))
    return added


def increase_model_bounds(model):
    current = 1000
    new_max = 99999999
    for reaction in model.reactions:
        if reaction.upper_bound == current:
            reaction.upper_bound = new_max
        if reaction.lower_bound == -current:
            reaction.lower_bound = -new_max


def map_reactions_list(map_path):
    """Extract reaction ids from map for FVA optimization

    :param map_path: string
    :return: list of strings
    """
    if not os.path.isfile(map_path):
        return []
    with open(map_path) as f:
        return [i['bigg_id'] for i in json.load(f)[1]['reactions'].values()]


def all_maps_reactions_list(model_name):
    """Extracts reaction ids from all the maps for the given model without duplicates

    :param model_name: string
    :return: list of strings
    """
    result = set()
    dirpath = MAPS_DIR + '/' + model_name
    for path, _, files in os.walk(dirpath):
        result.update(set([i for f in files for i in map_reactions_list(dirpath + '/' + f)]))
    return list(result)


class Response(object):
    def __init__(self, model, message, cache=None):
        self.model = model
        self.message = message
        self.method_name = message.get(SIMULATION_METHOD, 'fba')
        self.cache = cache
        if self.method_name in {'fva', 'pfba-fva'}:
            try:
                solution = self.solve_fva()
            except OptimizationError:
                logger.info('infeasible model for fva')
                self.flux = {}
                self.growth = 0.0
            else:
                df = solution.rename(index=str, columns={'maximum': 'upper_bound', 'minimum': 'lower_bound'})
                for key in ['lower_bound', 'upper_bound']:
                    df[key] = df[key].astype('float')
                self.flux = df.T.to_dict()
                self.growth = self.flux[MODEL_GROWTH_RATE[model.id]]['upper_bound']
        else:
            try:
                solution = self.solve()
            except OptimizationError:
                logger.info('infeasible model, returning measured fluxes only')
                changes = model.notes['changes']
                self.flux = {}
                self.growth = 0.0
                if 'measured' in changes:
                    ids_measured_reactions = set(rxn['id'] for rxn in changes['measured']['reactions'])
                    self.flux = dict((rxn.id, (rxn.upper_bound + rxn.lower_bound) / 2)
                                     for rxn in model.reactions if rxn.id in ids_measured_reactions)
            else:
                self.flux = solution.fluxes.to_dict()
                self.growth = self.flux[MODEL_GROWTH_RATE[model.id]]

    def solve_fva(self):
        fva_reactions = None
        if MAP in self.message:
            reaction_ids = map_reactions_list('{0}/{1}/{1}.{2}.json'.format(
                MAPS_DIR,
                self.model.id,
                self.message[MAP]
            ))
            if reaction_ids:
                reactions = [i for i in reaction_ids
                             if self.model.reactions.has_id(i)]
                fva_reactions = list(set(
                    reactions + [MODEL_GROWTH_RATE[self.model.id]]
                ))
        return METHODS[self.method_name](
            self.model,
            reactions=fva_reactions
        )

    def solve(self):
        t = time.time()
        if self.method_name in {'lmoma', 'room', 'moma'}:
            pfba_solution = pfba(self.model).fluxes.to_dict()
            if self.method_name == 'room':
                increase_model_bounds(self.model)
            solution = METHODS[self.method_name](self.model, cache=self.cache, reference=pfba_solution)
        else:
            solution = METHODS[self.method_name](self.model)
        logger.info('Model solved with method {} in {} sec'.format(self.method_name, time.time() - t))
        return solution

    def model_json(self):
        return model_to_dict(self.model)

    def fluxes(self):
        return self.flux

    async def theoretical_maximum_yield(self):
        res = {key: await phase_plane_to_dict(self.model, key) for key in self.message[OBJECTIVES]}
        logger.info(res)
        return res

    def growth_rate(self):
        return self.growth

    def measured_missing_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', deepcopy(EMPTY_CHANGES))['measured-missing']['reactions']]
        ))

    def removed_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', deepcopy(EMPTY_CHANGES))['removed']['reactions']]
        ))

    def added_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', deepcopy(EMPTY_CHANGES))['added']['reactions']
             if not is_dummy(i['id'])]
        ))


REQUEST_KEYS = [GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS]

APPLY_FUNCTIONS = {
    GENOTYPE_CHANGES: apply_genotype_changes,
    MEDIUM: apply_medium_changes,
    MEASUREMENTS: apply_measurement_changes,
}

RETURN_FUNCTIONS = {
    FLUXES: 'fluxes',
    TMY: 'theoretical_maximum_yield',
    MODEL: 'model_json',
    GROWTH_RATE: 'growth_rate',
    REMOVED_REACTIONS: 'removed_reactions',
    ADDED_REACTIONS: 'added_reactions',
    MISSING_MEASURED_REACTIONS: 'measured_missing_reactions',
}


async def modify_model(message, model):
    for key in REQUEST_KEYS:
        data = message.get(key, [])
        if data:
            modifications = await APPLY_FUNCTIONS[key](model, data)
            model = collect_changes(modifications)
    if REACTIONS_ADD in message:
        model = await apply_reactions_add(model, message[REACTIONS_ADD])
    if REACTIONS_KNOCKOUT in message:
        model = await apply_reactions_knockouts(model, message[REACTIONS_KNOCKOUT])
    return model


def collect_changes(modifications):
    model = modifications.model
    changes = model.notes.get('changes', deepcopy(EMPTY_CHANGES))
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


def restore_changes(model, changes):
    logger.info('Changes to restore: {}'.format(changes))
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


async def respond(message, model, db_key=None, cache=None):
    result = {}
    t = time.time()
    response = Response(model, message, cache=cache)
    for key in message['to-return']:
        if key == TMY:
            result[key] = await response.theoretical_maximum_yield()
        else:
            result[key] = getattr(response, RETURN_FUNCTIONS[key])()
    if db_key:
        result['model-id'] = db_key
    if REQUEST_ID in message:
        result[REQUEST_ID] = message[REQUEST_ID]
    logger.info('Response for {} is ready in {} sec'.format(message, time.time() - t))
    return result


async def model_ws_handler(request):
    ws = web.WebSocketResponse()
    model_id = request.match_info['model_id']
    cached_model = await restore_model(model_id)
    logger.info(model_from_changes.cache_info())
    if not cached_model:
        raise KeyError('No such model: {}'.format(model_id))
    model = cached_model.copy()
    model.notes = deepcopy(model.notes)
    cache = ProblemCache(model)
    await ws.prepare(request)
    try:
        async for msg in ws:
            logger.debug(msg)
            if msg.type == WSMsgType.TEXT:
                if msg.data == 'close':
                    await ws.close()
                else:
                    message = msg.json()
                    model = await modify_model(message, model)
                    ws.send_json(await respond(message, model, cache=cache))
            elif msg.type == WSMsgType.ERROR:
                logger.error('Websocket for model_id {} closed with exception {}'.format(model_id, ws.exception()))
    except asyncio.CancelledError:
        logger.info('Websocket for model_id {} cancelled'.format(model_id))
    await ws.close()
    return ws


async def model_handler(request):
    wild_type_id = request.match_info['model_id']
    data = await request.json()
    if 'message' not in data:
        return web.HTTPBadRequest()
    message = data['message']
    db_key = key_from_model_info(wild_type_id, message)
    model = await restore_from_db(db_key)
    logger.info(model_from_changes.cache_info())
    if not model:
        model = find_in_memory(wild_type_id)
        if not model:
            return web.HTTPNotFound()
        model = model.copy()
        model.notes = deepcopy(model.notes)
        model = await modify_model(message, model)
        db_key = await save_changes_to_db(model, wild_type_id, message)
    if message.get(SIMULATION_METHOD) == 'room':
        model = model.copy()
    return web.json_response(await respond(message, model, db_key))


async def maps_handler(request):
    return web.json_response(MAP_DICTIONARY)


async def model_options_handler(request):
    return web.json_response(SPECIES_TO_MODEL[request.match_info['species']])


async def map_handler(request):
    # FIXME: A malicious user can access any JSON file in the system this way.
    filepath = '{}/{}/{}.{}.json'.format(
        MAPS_DIR, request.GET['model'], request.GET['model'], request.GET['map']
    )
    with open(filepath) as f:
        return web.json_response(json.load(f))


app = web.Application()
app.router.add_route('GET', '/wsmodels/{model_id}', model_ws_handler)
app.router.add_route('GET', '/maps', maps_handler)
app.router.add_route('GET', '/map', map_handler)
app.router.add_route('GET', '/model-options/{species}', model_options_handler)
app.router.add_route('POST', '/models/{model_id}', model_handler)

# Configure default CORS settings.
cors = aiohttp_cors.setup(app, defaults={
    "*": aiohttp_cors.ResourceOptions(
        allow_credentials=True,
        expose_headers="*",
        allow_headers="*",
    )
})

# Configure CORS on all routes.
for route in list(app.router.routes()):
    cors.add(route)


async def start(loop):
    await loop.create_server(app.make_handler(), '0.0.0.0', 8000)
    logger.info('Web server is up')


if __name__ == '__main__':
    loop = asyncio.get_event_loop()
    loop.run_until_complete(start(loop))
    try:
        loop.run_forever()
    except KeyboardInterrupt:
        pass
