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

import logging

from cobra.io.dict import model_to_dict
from flask import Response, jsonify, request
from prometheus_client import CONTENT_TYPE_LATEST, CollectorRegistry, generate_latest
from prometheus_client.multiprocess import MultiProcessCollector

import model.constants as constants
from aiohttp import WSMsgType, web
from model.operations import modify_model
from model.response import respond
from model.storage import Models, key_from_model_info, model_from_changes, restore_from_db, save_changes_to_db


logger = logging.getLogger(__name__)


def model(model_id):
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    wild_type_id = model_id
    if 'message' not in request.json:
        return "Missing field 'message'", 400
    message = request.json['message']

    mutated_model_id = key_from_model_info(wild_type_id, message)
    model = restore_from_db(mutated_model_id)
    logger.info(model_from_changes.cache_info())
    if not model:
        model = Models.get(wild_type_id)
        if not model:
            return f"Model '{wild_type_id}' not found", 404
        model = model.copy()
        model = modify_model(message, model)
        mutated_model_id = save_changes_to_db(model, wild_type_id, message)
    return jsonify(respond(model, message, mutated_model_id=mutated_model_id))


def model_get(model_id):
    wild_type_id = model_id
    wild_type_model = Models.get(wild_type_id)
    return jsonify(model_to_dict(wild_type_model))


def model_diff(model_id):
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    wild_type_id = model_id
    if 'message' not in request.json:
        return "Missing field 'message'", 400
    message = request.json['message']

    wild_type_model = Models.get(wild_type_id)
    if not wild_type_model:
        return f"Model '{wild_type_id}' not found", 404

    mutated_model_id = key_from_model_info(wild_type_id, message, version=1)
    mutated_model = restore_from_db(mutated_model_id)
    logger.info(model_from_changes.cache_info())

    if not mutated_model:
        # @matyasfodor Not sure about how this will perform, copy is an expensive
        # operation. We could just modify the previous state (get it from LRU cache)
        # but then we'd modify the returned object, it should not be accessed more
        # than once. One idea is to delete the retrieved obejct (only use it once).
        # It is possible with https://pypi.python.org/pypi/pylru
        mutated_model = wild_type_model.copy()
        mutated_model = modify_model(message, mutated_model)
        mutated_model_id = save_changes_to_db(mutated_model, wild_type_id, message, version=1)

    diff = respond(mutated_model, message, mutated_model_id, wild_type_model_id=wild_type_id)
    return jsonify(diff)


def model_info(model_id):
    try:
        wild_type_model = Models.get(model_id)
        medium = [{
            'id': reaction_id,
            'name': wild_type_model.reactions.get_by_id(reaction_id).name.replace('exchange', '').strip()
        } for reaction_id in wild_type_model.medium]
        return jsonify({'medium': medium})
    except KeyError:
        return f"Unknown model {model_id}", 400


def model_options(species):
    try:
        return jsonify(constants.SPECIES_TO_MODEL[species])
    except KeyError:
        return f"Unknown species {species}", 400


def metrics():
    return Response(generate_latest(MultiProcessCollector(CollectorRegistry())), mimetype=CONTENT_TYPE_LATEST)
