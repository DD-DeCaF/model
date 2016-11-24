import asyncio
import aiohttp_cors
import gnomic
import json
import hashlib
import time
import aioredis
import os
from copy import deepcopy
from collections import namedtuple
from functools import lru_cache
from aiohttp import web, WSMsgType
from cameo.data import metanetx
from cameo import load_model
from cameo import phenotypic_phase_plane
from cobra.io.json import to_json, reaction_to_dict, reaction_from_dict, gene_to_dict, \
    metabolite_to_dict, metabolite_from_dict
from driven.generic.adapter import get_existing_metabolite, GenotypeChangeModel, MediumChangeModel, \
    MeasurementChangeModel, full_genotype, feature_id
from venom.rpc.comms.grpc import Client
from model import logger
from model.messages import GeneToReactionsRemote, GeneRequest


ORGANISMS = {
    'iJO1366',
    'iMM904',
    'iMM1415',
    'iNJ661',
    'e_coli_core',
}


async def redis_client():
    return await aioredis.create_redis((os.environ['REDIS_PORT_6379_TCP_ADDR'], 6379), loop=asyncio.get_event_loop())


class Models(object):
    MODELS = {
        v: load_model(v) for v in ORGANISMS
    }
    print('Models are ready')

GENOTYPE_CHANGES = 'genotype-changes'
MEDIUM = 'medium'
MEASUREMENTS = 'measurements'
REACTIONS = 'reactions-knockout'
REACTIONS_UNDO = 'reactions-knockout-undo'
MODEL = 'model'
FLUXES = 'fluxes'
GROWTH_RATE = 'growth-rate'
TMY = 'tmy'
OBJECTIVES = 'objectives'
REQUEST_ID = 'request-id'
REMOVED_REACTIONS = 'removed-reactions'

EMPTY_CHANGES = {
    'added': {
        'reactions': [],
        'metabolites': [],
    },
    'removed': {
        'genes': [],
        'reactions': [],
    }
}

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


@lru_cache(maxsize=16)
def model_from_changes(changes):
    changes = json.loads(changes)
    model = find_in_memory(changes['model']).copy()
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


def key_from_model_info(model_id, message):
    """Generate hash string from model information which will later be used as db key

    :param model_id: str
    :param message: dict
    :return: str
    """
    d = {k: message.get(k, []) for k in {GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS}}
    d['model_id'] = model_id
    return hashlib.sha224(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()


async def save_changes_to_db(model, model_id, message):
    """Store model in cache database

    :param model: Cameo model
    :param model_id: str
    :param message: dict
    :return: str (cache database key)
    """
    db_key = key_from_model_info(model_id, message)
    redis = await redis_client()
    await redis.set(db_key, json.dumps({'model': model.id, 'changes': model.notes.get('changes', deepcopy(EMPTY_CHANGES))}))
    redis.close()
    logger.info('Model created on the base of {} with message {} saved as {}'.format(model_id, message, db_key))
    return db_key


class NoIDMapping(Exception):
    def __init__(self, metabolite_id):
        self.value = metabolite_id

    def __str__(self):
        return 'No Metanetx mapping for metabolite {}'.format(self.value)


def existing_metabolite(model, metabolite_id):
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
    return get_existing_metabolite(
        mnx_id,
        model,
        '_e'
    )


def product_reaction_variable(model, metabolite_id):
    """Find a medium exchange reaction in the model for the given metabolite id

    :param model: cameo model
    :param metabolite_id: string of format <database>:<id>, f.e. chebi:12345
    :return:
    """
    metabolite = existing_metabolite(model, metabolite_id)
    if not metabolite:
        return None
    exchange_reactions = list(set(metabolite.reactions).intersection(model.exchanges))
    if len(exchange_reactions) != 1:
        logger.info('More than one exchange reaction in model {} for metabolite {}'.format(model.id, metabolite_id))
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
    result = phenotypic_phase_plane(model, [reaction]).data_frame.to_dict()
    for k, v in result.items():
        result[k] = [v[point] for point in sorted(v.keys())]
    return result


async def apply_genotype_changes(model, genotype_changes):
    """Apply genotype changes to cameo model.

    :param initial_model: cameo model
    :param genotype_changes: list of strings, f.e. ['-tyrA::kanMX+', 'kanMX-']
    :return:
    """
    logger.info('Genotype changes {}'.format(genotype_changes))
    genotype_features = full_genotype(genotype_changes)
    genes_to_reactions = await call_genes_to_reactions(genotype_features)
    logger.info('Genes to reaction: {}'.format(genes_to_reactions))
    return GenotypeChangeModel(model, genotype_features, genes_to_reactions)


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
    client = Client(GeneToReactionsRemote, '139.59.133.210', 50053)
    identifiers = list(new_features_identifiers(genotype_features))
    results = await asyncio.gather(*[
        client.invoke(
            GeneToReactionsRemote.reactions,
            GeneRequest(gene=identifier)
        ) for identifier in identifiers
        ])
    return {k: dict(zip(*(v.reactions_ids, v.equations))) for k, v in zip(identifiers, results)}


def convert_mg_to_mmol(mg, formula_weight):
    return mg * (1/formula_weight)


async def apply_measurement_changes(model, measurements):
    measurements = convert_measurements_to_mmol(measurements, model)
    return MeasurementChangeModel(model, measurements)


def convert_measurements_to_mmol(measurements, model):
    for i in range(len(measurements)):
        if 'unit' not in measurements[i]:
            continue
        if measurements[i]['unit'] == 'mg':
            metabolite = existing_metabolite(model, measurements[i]['id'])
            if metabolite:
                measurements[i]['measurement'] = convert_mg_to_mmol(
                    measurements[i]['measurement'],
                    metabolite.formula_weight
                )
                measurements[i]['unit'] = 'mmol'
                logger.info('Converted metabolite {} from mg to mmol'.format(measurements[i]['id']))
    return measurements


async def apply_medium_changes(model, medium):
    return MediumChangeModel(model, medium)


ReactionKnockouts = namedtuple('ReactionKnockouts', ['model', 'changes'])

async def apply_reactions_knockouts(model, reactions_ids):
    reactions = [model.reactions.get_by_id(r_id).copy() for r_id in reactions_ids]
    for r in reactions:
        model.reactions.get_by_id(r.id).knock_out()
    return ReactionKnockouts(model, {'removed': {'reactions': reactions}})


async def undo_reactions_knockouts(model, reaction_ids):
    reaction_ids = set(reaction_ids)
    logger.info('Reactions to undo knockout: {}'.format(reaction_ids))
    if not model.notes.get('changes'):
        return ReactionKnockouts(model, {})
    reactions = model.notes['changes']['removed']['reactions']
    for reaction in reactions:
        if reaction['id'] in reaction_ids:
            model.reactions.get_by_id(reaction['id']).change_bounds(
                lb=reaction['lower_bound'],
                ub=reaction['upper_bound']
            )
    model.notes['changes']['removed']['reactions'] = [i for i in reactions if i['id'] not in reaction_ids]
    return ReactionKnockouts(model, {})


class Response(object):
    def __init__(self, model, message):
        self.model = model
        self.message = message
        solution = self.model.solve()
        self.flux = solution.fluxes
        self.growth = solution.objective_value

    def model_json(self):
        return json.loads(to_json(self.model))

    def fluxes(self):
        return self.flux

    def theoretical_maximum_yield(self):
        return {key: phase_plane_to_dict(self.model, key) for key in self.message[OBJECTIVES]}

    def growth_rate(self):
        return self.growth

    def removed_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', deepcopy(EMPTY_CHANGES))['removed']['reactions']]
        ))


REQUEST_KEYS = [GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS, REACTIONS, REACTIONS_UNDO]


APPLY_FUNCTIONS = {
    GENOTYPE_CHANGES: apply_genotype_changes,
    MEDIUM: apply_medium_changes,
    MEASUREMENTS: apply_measurement_changes,
    REACTIONS: apply_reactions_knockouts,
    REACTIONS_UNDO: undo_reactions_knockouts,
}

RETURN_FUNCTIONS = {
    FLUXES: 'fluxes',
    TMY: 'theoretical_maximum_yield',
    MODEL: 'model_json',
    GROWTH_RATE: 'growth_rate',
    REMOVED_REACTIONS: 'removed_reactions',
}

async def modify_model(message, model):
    for key in REQUEST_KEYS:
        data = message.get(key, [])
        if data:
            modifications = await APPLY_FUNCTIONS[key](model, data)
            model = collect_changes(modifications)
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
            changes[action][entity].extend([to_dict[entity](i) for i in value])
    model.notes['changes'] = changes
    return model


def restore_changes(model, changes):
    model = apply_additions(model, changes['added'])
    model = apply_removals(model, changes['removed'])
    return model


def apply_additions(model, changes):
    model = add_metabolites(model, changes['metabolites'])
    model = add_reactions(model, changes['reactions'])
    return model


def add_reactions(model, changes):
    reactions = [reaction_from_dict(r, model) for r in changes]
    to_add = []
    for reaction in reactions:
        if model.reactions.has_id(reaction.id):
            model.reactions.get_by_id(reaction.id).change_bounds(lb=reaction.lower_bound, ub=reaction.upper_bound)
        else:
            to_add.append(reaction)
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


def respond(message, model, db_key=None):
    result = {}
    response = Response(model, message)
    for key in message['to-return']:
        result[key] = getattr(response, RETURN_FUNCTIONS[key])()
    if db_key:
        result['model-id'] = db_key
    if REQUEST_ID in message:
        result[REQUEST_ID] = message[REQUEST_ID]
    logger.info('Response for {} is ready'.format(message))
    return result


async def model_ws_handler(request):
    ws = web.WebSocketResponse()
    model_id = request.match_info['model_id']
    model = await restore_model(model_id)
    logger.info(model_from_changes.cache_info())
    if not model:
        raise KeyError('No such model: {}'.format(model_id))
    model = model.copy()
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
                    ws.send_json(respond(message, model))
            elif msg.type == WSMsgType.ERROR:
                logger.error('Websocket for model_id {} closed with exception {}'.format(model_id, ws.exception()))
    except asyncio.CancelledError:
        logger.info('Websocket for model_id {} cancelled'.format(model_id))
    await ws.close()
    return ws


async def model_handler(request):
    model_id = request.match_info['model_id']
    data = await request.json()
    if 'message' not in data:
        return web.HTTPBadRequest()
    message = data['message']
    db_key = key_from_model_info(model_id, message)
    model = await restore_from_db(db_key)
    logger.info(model_from_changes.cache_info())
    if not model:
        model = find_in_memory(model_id)
        if not model:
            return web.HTTPNotFound()
        model = await modify_model(message, model.copy())
        db_key = await save_changes_to_db(model, model_id, message)
    return web.json_response(respond(message, model, db_key))


app = web.Application()
app.router.add_route('GET', '/wsmodels/{model_id}', model_ws_handler)
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
