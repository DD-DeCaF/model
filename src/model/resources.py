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

from model import constants, storage
from model.adapter import adapt_from_medium, adapt_from_genotype, adapt_from_measurements
from model.operations import modify_model
from model.simulations import simulate


logger = logging.getLogger(__name__)


def species(species):
    return jsonify([model_meta.model_id for model_meta in storage.MODELS if model_meta.species == species])


def delta_create():
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    try:
        model_id = request.json['model_id']
    except KeyError:
        return "Missing field 'model_id'", 400

    try:
        model_meta = storage.get(model_id)
    except KeyError:
        return f"Unknown model '{model_id}'", 400

    # Build list of operations to perform on the model
    operations = []
    if 'medium' in request.json:
        operations.extend(adapt_from_medium(model_meta.model, request.json['medium']))

    if 'genotype' in request.json:
        operations.extend(adapt_from_genotype(model_meta.model, request.json['genotype']))

    if 'measurements' in request.json:
        operations.extend(adapt_from_measurements(model_meta.model, request.json['measurements']))

    if 'operations' in request.json:
        operations.extend(request.json['operations'])

    # TODO: Save operations, get id
    def create_delta(operations):
        pass

    operations_id = create_delta(operations)
    return jsonify({'id': operations_id, 'operations': operations})


def delta_get(delta_id):
    pass


def delta_update(delta_id):
    pass


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

    # Parse solver request args and set defaults
    method = message.get(constants.SIMULATION_METHOD, 'fba')
    objective_id = message.get(constants.OBJECTIVE)
    objective_direction = message.get(constants.OBJECTIVE_DIRECTION)
    tmy_objectives = message.get(constants.TMY_OBJECTIVES, [])
    to_return = message.get('to-return')

    result = simulate(model, method, objective_id, objective_direction, tmy_objectives, to_return)
    response = {key: value for key, value in result.items() if key in to_return}
    response['model-id'] = changes_key
    return jsonify(response)


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


def healthz():
    """
    HTTP endpoint for readiness probes.

    Return an empty response. This response will not be ready until the application has finished initializing, e.g.,
    preloading models, which takes a few minutes.
    """
    return ""
