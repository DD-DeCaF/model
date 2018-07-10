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

from model import storage
from model.operations import modify_model
from model.response import respond


logger = logging.getLogger(__name__)


def species(species):
    try:
        return jsonify([model_meta.model_id for model_meta in storage.MODELS if model_meta.species == species])
    except KeyError:
        return f"Unknown species {species}", 404


def model_get(model_id):
    try:
        return jsonify(model_to_dict(storage.get(model_id).model))
    except KeyError:
        return f"Unknown model {model_id}", 404


def model_modify_simulate(model_id):
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    if 'message' not in request.json:
        return "Missing field 'message'", 400
    message = request.json['message']

    try:
        model = storage.restore_from_message(model_id, message)
        changes_key = storage._changes_key(model_id, message)  # FIXME: shouldn't need to know changes_key here
    except KeyError:
        try:
            model_meta = storage.get(model_id)
        except KeyError:
            return f"Unknown model {model_id}", 404
        model = model_meta.model.copy()
        model = modify_model(message, model)
        changes_key = storage.save_changes(model, message)
    return jsonify(respond(model, message, mutated_model_id=changes_key))


def model_medium(model_id):
    try:
        model = storage.get(model_id).model
        medium = [{
            'id': reaction_id,
            'name': model.reactions.get_by_id(reaction_id).name.replace('exchange', '').strip()
        } for reaction_id in model.medium]
        return jsonify({'medium': medium})
    except KeyError:
        return f"Unknown model {model_id}", 404


def metrics():
    return Response(generate_latest(MultiProcessCollector(CollectorRegistry())), mimetype=CONTENT_TYPE_LATEST)
