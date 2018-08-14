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

from redis import Redis

from model.app import app


logger = logging.getLogger(__name__)
redis = Redis(app.config['REDIS_ADDR'], app.config['REDIS_PORT'])


def save(model_id, conditions, operations):
    """Save adapted changes in cache based on the given model and list of operations"""
    delta_id = _delta_key(model_id, conditions)
    redis.set(delta_id, json.dumps(operations))
    logger.info(f"Stored changes for {model_id} as '{delta_id}'")
    return delta_id


def load(model_id, conditions):
    """
    Return the pre-computed, cached operations based on the given conditions, if it exists.

    :param model_id: cobra model id
    :param conditions: the raw conditions and omics applied to yield operations
    :raises: KeyError if the given input does not have a cached result
    :return: cached operations
    """
    return load_from_key(_delta_key(model_id, conditions))


def load_from_key(delta_id):
    """Restore model with modifications from cache based on the given model and message"""
    operations_json = redis.get(delta_id)
    if operations_json is None:
        logger.debug(f"Changes '{delta_id}' not found in cache")
        raise KeyError(f"Changes '{delta_id}' not found in cache")
    return json.loads(operations_json.decode('utf-8'))


def _delta_key(model_id, conditions):
    """Generate unique input string used for cache key based on the model id and given conditions"""
    unique_key_input = json.dumps({'model_id': model_id, 'conditions': conditions}, sort_keys=True)
    return hashlib.sha224(unique_key_input.encode('utf-8')).hexdigest()
