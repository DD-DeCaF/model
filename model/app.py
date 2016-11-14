import asyncio
import aiohttp_cors
import gnomic
from aiohttp import web, WSMsgType
from cameo.data import metanetx
from cameo import load_model
from cameo import phenotypic_phase_plane
from driven.generic.adapter import get_existing_metabolite, GenotypeChangeModel, MediumChangeModel, \
    MeasurementChangeModel, full_genotype, feature_id
from venom.rpc.comms.grpc import Client
from model import logger
from model.messages import GeneToReactionsRemote, GeneRequest

ORGANISMS = {
    'ECO': 'iJO1366',
    'SCE': 'iMM904',
    'CHO': 'iMM1415',
    'COG': 'iNJ661',
}


class Models(object):
    MODELS = {
        k: load_model(v) for k, v in ORGANISMS.items()
    }
    print('Models are ready')

GENOTYPE_CHANGES = 'genotype-changes'
MEDIUM = 'medium'
MEASUREMENTS = 'measurements'
REACTIONS = 'reactions-knockout'
FLUXES = 'fluxes'
TMY = 'tmy'
OBJECTIVES = 'objectives'


def restore_model(model_id):
    return Models.MODELS.get(model_id).copy()


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
    return GenotypeChangeModel(model, genotype_features, genes_to_reactions).model


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


async def apply_measurement_changes(model, measurements):
    return MeasurementChangeModel(model, measurements).model


async def apply_medium_changes(model, medium):
    return MediumChangeModel(model, medium).model


async def apply_reactions_knockouts(model, reactions_ids):
    for r_id in reactions_ids:
        reaction = model.reactions.get_by_id(r_id)
        reaction.knock_out()
    return model


def fluxes(model, message):
    return model.solve().fluxes


def theoretical_maximum_yield(model, message):
    return {key: phase_plane_to_dict(model, key) for key in message[OBJECTIVES]}


REQUEST_KEYS = [GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS, REACTIONS]


APPLY_FUNCTIONS = {
    GENOTYPE_CHANGES: apply_genotype_changes,
    MEDIUM: apply_medium_changes,
    MEASUREMENTS: apply_measurement_changes,
    REACTIONS: apply_reactions_knockouts,
}

RETURN_FUNCTIONS = {
    FLUXES: fluxes,
    TMY: theoretical_maximum_yield,
}

async def modify_model(message, model):
    for key in REQUEST_KEYS:
        data = message.get(key, [])
        if data:
            model = await APPLY_FUNCTIONS[key](model, data)
    return model


def respond(message, model):
    result = {}
    for key in message['to-return']:
        result[key] = RETURN_FUNCTIONS[key](model, message)
    return result


async def model_handler(request):
    ws = web.WebSocketResponse()
    model_id = request.match_info['model_id']
    model = restore_model(model_id)
    if not model:
        raise KeyError('No such model: {}'.format(model_id))
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


app = web.Application()
app.router.add_route('GET', '/models/{model_id}', model_handler)


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
