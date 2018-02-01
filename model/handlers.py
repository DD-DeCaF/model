import asyncio
from aiohttp import web, WSMsgType
import json
import logging
import os

from cobra.io.dict import model_to_dict

import model.constants as constants
from model.storage import (restore_model, model_from_changes, key_from_model_info,
                           restore_from_db, Models, save_changes_to_db)
from model.operations import modify_model
from model.response import respond

LOGGER = logging.getLogger(__name__)

async def model_ws_full(request):
    return await model_ws(request, False)

async def model_ws_json_diff(request):
    return await model_ws(request, True)

async def model_ws(request, diff=False):
    ws = web.WebSocketResponse()
    model_id = request.match_info['model_id']
    cached_model = await restore_model(model_id)
    LOGGER.info(model_from_changes.cache_info())
    if not cached_model:
        raise KeyError('No such model: {}'.format(model_id))
    model = cached_model.copy()
    await ws.prepare(request)
    try:
        async for msg in ws:
            LOGGER.debug(msg)
            if msg.type == WSMsgType.TEXT:
                if msg.data == 'close':
                    await ws.close()
                else:
                    message = msg.json()
                    model = await modify_model(message, model)
                    await ws.send_json(await respond(model, message, wild_type_model_id=model_id if diff else None))
            elif msg.type == WSMsgType.ERROR:
                LOGGER.error('Websocket for model_id %s closed with exception %s', model_id, ws.exception())
    except asyncio.CancelledError as ex:
        LOGGER.debug('Websocket for model_id %s cancelled, %s', model_id, str(ex))
    await ws.close()
    return ws


async def model(request):
    wild_type_id = request.match_info['model_id']
    data = await request.json()

    try:
        message = data['message']
    except KeyError:
        return web.HTTPBadRequest()

    mutated_model_id = key_from_model_info(wild_type_id, message)
    model = await restore_from_db(mutated_model_id)
    LOGGER.info(model_from_changes.cache_info())
    if not model:
        model = Models.get(wild_type_id)
        if not model:
            return web.HTTPNotFound()
        model = model.copy()
        model = await modify_model(message, model)
        mutated_model_id = await save_changes_to_db(model, wild_type_id, message)
    return web.json_response(await respond(model, message, mutated_model_id=mutated_model_id))


async def model_get(request):
    wild_type_id = request.match_info['model_id']
    wild_type_model = Models.get(wild_type_id)
    return web.json_response(model_to_dict(wild_type_model))


async def model_diff(request):
    wild_type_id = request.match_info['model_id']
    data = await request.json()

    try:
        message = data['message']
    except KeyError:
        return web.HTTPBadRequest()

    wild_type_model = Models.get(wild_type_id)
    if not wild_type_model:
        return web.HTTPNotFound()
    mutated_model_id = key_from_model_info(wild_type_id, message, version=1)
    mutated_model = await restore_from_db(mutated_model_id)
    LOGGER.info(model_from_changes.cache_info())

    if not mutated_model:
        # @matyasfodor Not sure about how this will perform, copy is an expensive
        # operation. We could just modify the previous state (get it from LRU cache)
        # but then we'd modify the returned object, it should not be accessed more
        # than once. One idea is to delete the retrieved obejct (only use it once).
        # It is possible with https://pypi.python.org/pypi/pylru
        mutated_model = wild_type_model.copy()
        mutated_model = await modify_model(message, mutated_model)
        mutated_model_id = await save_changes_to_db(mutated_model, wild_type_id, message, version=1)

    diff = await respond(mutated_model, message, mutated_model_id, wild_type_model_id=wild_type_id)
    return web.json_response(diff)


async def model_info(request):
    wild_type_id = request.match_info['model_id']
    wild_type_model = Models.get(wild_type_id)
    medium = [{
        'id': reaction_id,
        'name': wild_type_model.reactions.get_by_id(reaction_id).name.replace('exchange', '').strip()
    } for reaction_id in wild_type_model.medium]
    return web.json_response({'medium': medium})


async def maps(request):
    return web.json_response(constants.MAP_DICTIONARY)


async def model_options(request):
    return web.json_response(constants.SPECIES_TO_MODEL[request.match_info['species']])


# TODO @matyasfodor This really shouldn't live here. This is basically a static file service
async def map(request):
    modelId = request.GET['model']
    mapId = request.GET['map']
    directory = os.path.realpath(constants.MAPS_DIR)
    filepath = os.path.join(directory, modelId, '{}.{}.json'.format(modelId, mapId))
    if os.path.commonprefix((os.path.realpath(filepath), directory)) != directory:
        return web.HTTPBadRequest()
    try:
        with open(filepath) as f:
            return web.json_response(json.load(f))
    except FileNotFoundError:
        return web.HTTPNotFound()
