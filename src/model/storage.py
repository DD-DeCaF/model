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

import hashlib
import json
import logging
import os
from functools import lru_cache

from cobra.io import model_to_dict, read_sbml_model
from redis import Redis

from model import constants
from model.app import app
from model.operations import restore_changes
from model.utils import log_time


logger = logging.getLogger(__name__)
redis = Redis(app.config['REDIS_ADDR'], app.config['REDIS_PORT'])


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


def save_changes_to_db(model, wild_type_id, message, version=None):
    """Store model in cache database

    :param model: Cameo model
    :param wild_type_id: str
    :param message: dict
    :return: mutated_model_id (cache database key)
    """
    mutated_model_id = key_from_model_info(wild_type_id, message, version)

    value = json.dumps({'model': model.id, 'changes': model.notes.get('changes', constants.get_empty_changes())})
    redis.set(mutated_model_id, value)
    logger.info(f"Model created on the base of {wild_type_id} with message {message} saved as {mutated_model_id}")
    return mutated_model_id


def read_model(model_id):
    with log_time(operation=f"Read model {model_id} from SBML file"):
        model = read_sbml_model(os.path.join('data', 'models', model_id + '.sbml.gz'))
        model.solver = 'cplex'
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
    logger.info("Preloading models")
    with log_time(operation="Preload models"):
        for model_id in constants.MODELS:
            Models.get_dict(model_id)


if app.config['ENVIRONMENT'] in ['production', 'staging']:
    preload_cache()


def find_changes_in_db(model_id):
    dumped_changes = redis.get(model_id)
    if dumped_changes is not None:
        return dumped_changes.decode('utf-8')
    return None


def restore_from_db(model_id):
    with log_time(operation=f"Restored model {model_id} from db"):
        changes = find_changes_in_db(model_id)
        if not changes:
            return None
        return model_from_changes(changes)


@lru_cache(maxsize=2 ** 6)
def model_from_changes(changes):
    changes = json.loads(changes)
    model = Models.get(changes['model']).copy()
    model = restore_changes(model, changes['changes'])
    model.notes['changes'] = changes['changes']
    return model


def restore_model(model_id):
    """Try to restore model by model id.
    NOTE: if model is found in memory, the original model is returned - to modify, make a copy

    :param model_id: str
    :return: Cameo model or None
    """
    model = Models.get(model_id)
    if model:
        logger.info('Wild type model with id %s is found', model_id)
        return model
    model = restore_from_db(model_id)
    if model:
        logger.info('Model with id %s found in database', model_id)
        return model
    logger.info('No model with id %s', model_id)
    return None
