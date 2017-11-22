import aioredis
import asyncio
from functools import lru_cache
import hashlib
import json
import logging
import os
import pickle
import time

from cobra.io import read_sbml_model, model_to_dict

import model.constants as constants
from model.utils import timing
from model.operations import restore_changes


LOGGER = logging.getLogger(__name__)


def key_from_model_info(wild_type_id, message, version=None):
    """Generate hash string from model information which will later be used as db key

    :param wild_type_id: str
    :param message: dict
    :return: str
    """
    d = {k: message.get(k, []) for k in constants.MESSAGE_HASH_KEYS}
    d['model_id'] = wild_type_id

    if version:
        d['version'] = version

    return hashlib.sha224(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()


async def redis_client():
    return await aioredis.create_redis((os.environ['REDIS_ADDR'], os.environ['REDIS_PORT']),
                                       loop=asyncio.get_event_loop())


async def save_changes_to_db(model, wild_type_id, message, version=None):
    """Store model in cache database

    :param model: Cameo model
    :param wild_type_id: str
    :param message: dict
    :return: mutated_model_id (cache database key)
    """
    mutated_model_id = key_from_model_info(wild_type_id, message, version)

    value = json.dumps({'model': model.id, 'changes': model.notes.get('changes', constants.get_empty_changes())})
    redis = await redis_client()
    with (await redis) as connection:
        await connection.set(mutated_model_id, value)
    LOGGER.info('Model created on the base of %s with message %s saved as %s', wild_type_id, message, mutated_model_id)
    return mutated_model_id


def read_model(model_id):
    model = read_sbml_model(os.path.join(os.path.dirname(__file__), 'data', model_id + '.sbml.gz'))
    model.notes['namespace'] = constants.MODEL_NAMESPACE[model_id]
    return model


class Models(object):
    @classmethod
    @lru_cache()
    def get(cls, model_id):
        try:
            return read_model(model_id)
        except FileNotFoundError:
            return None

    @classmethod
    @lru_cache()
    def get_dict(cls, model_id):
        model = cls.get(model_id)
        if model is not None:
            return model_to_dict(model)
        return None


def preload_cache():
    for model_id in constants.MODELS:
        Models.get_dict(model_id)


if constants.ENV == constants.ENV_PROD:
    preload_cache()


async def find_changes_in_db(model_id):
    dumped_changes = None
    redis = await redis_client()
    with (await redis) as connection:
        dumped_changes = await connection.get(model_id)

    if dumped_changes is not None:
        return dumped_changes.decode('utf-8')
    return None


async def restore_from_db(model_id):
    t = time.time()
    changes = await find_changes_in_db(model_id)
    if not changes:
        return None
    model = model_from_changes(changes)
    t = time.time() - t
    LOGGER.info('Model with db key %s is ready in %s sec', model_id, t)
    return model


@lru_cache(maxsize=2 ** 6)
def model_from_changes(changes):
    changes = json.loads(changes)
    model = Models.get(changes['model']).copy()
    model = restore_changes(model, changes['changes'])
    model.notes['changes'] = changes['changes']
    return model


async def restore_model(model_id):
    """Try to restore model by model id.
    NOTE: if model is found in memory, the original model is returned - to modify, make a copy

    :param model_id: str
    :return: Cameo model or None
    """
    model = Models.get(model_id)
    if model:
        LOGGER.info('Wild type model with id %s is found', model_id)
        return model
    model = await restore_from_db(model_id)
    if model:
        LOGGER.info('Model with id %s found in database', model_id)
        return model
    LOGGER.info('No model with id %s', model_id)
    return None
