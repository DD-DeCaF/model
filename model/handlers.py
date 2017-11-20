import asyncio
from aiohttp import web, WSMsgType
import json

from cobra.io.dict import model_to_dict

import model.consts as consts
from model.logger import logger
from model.memory_storage import restore_model, model_from_changes, key_from_model_info, restore_from_db, find_in_memory, save_changes_to_db
from model.model_operations import modify_model
from model.response import respond

async def model_ws(request):
    ws = web.WebSocketResponse()
    model_id = request.match_info['model_id']
    cached_model = await restore_model(model_id)
    logger.info(model_from_changes.cache_info())
    if not cached_model:
        raise KeyError('No such model: %s', model_id)
    model = cached_model.copy()
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
                    ws.send_json(await respond(model, message))
            elif msg.type == WSMsgType.ERROR:
                logger.error('Websocket for model_id %s closed with exception %s', model_id, ws.exception())
    except asyncio.CancelledError:
        logger.info('Websocket for model_id %s cancelled', model_id)
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
    logger.info(model_from_changes.cache_info())
    if not model:
        model = find_in_memory(wild_type_id)
        if not model:
            return web.HTTPNotFound()
        model = model.copy()
        model = await modify_model(message, model)
        mutated_model_id = await save_changes_to_db(model, wild_type_id, message)
    return web.json_response(await respond(model, message, mutated_model_id))


async def model_get(request):
    wild_type_id = request.match_info['model_id']
    wild_type_model = find_in_memory(wild_type_id)
    return web.json_response(model_to_dict(wild_type_model))


async def model_diff(request):
    wild_type_id = request.match_info['model_id']
    data = await request.json()

    try:
        message = data['message']
    except KeyError:
        return web.HTTPBadRequest()

    wild_type_model = find_in_memory(wild_type_id)
    if not wild_type_model:
        return web.HTTPNotFound()
    mutated_model_id = key_from_model_info(wild_type_id, message, version=1)
    mutated_model = await restore_from_db(mutated_model_id)
    logger.info(model_from_changes.cache_info())

    if not mutated_model:
        # @matyasfodor Not sure about how this will perform, copy is an expensive
        # operation. We could just modify the previous state (get it from LRU cache)
        # but then we'd modify the returned object, it should not be accessed more
        # than once. One idea is to delete the retrieved obejct (only use it once).
        # It is possible with https://pypi.python.org/pypi/pylru
        mutated_model = wild_type_model.copy()
        mutated_model = await modify_model(message, mutated_model)
        mutated_model_id = await save_changes_to_db(mutated_model, wild_type_id, message, version=1)

    diff = await respond(mutated_model, message, mutated_model_id, wild_type_model)
    return web.json_response(diff)

async def maps(request):
    return web.json_response(consts.MAP_DICTIONARY)


async def model_options(request):
    return web.json_response(consts.SPECIES_TO_MODEL[request.match_info['species']])


async def map(request):
    # FIXME: A malicious user can access any JSON file in the system this way.
    filepath = '{}/{}/{}.{}.json'.format(
        consts.MAPS_DIR, request.GET['model'], request.GET['model'], request.GET['map']
    )
    with open(filepath) as f:
        return web.json_response(json.load(f))
