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

from model import deltas, storage
from model.adapter import adapt_from_genotype, adapt_from_measurements, adapt_from_medium
from model.operations import apply_operations
from model.simulations import simulate


logger = logging.getLogger(__name__)


def species(species):
    return jsonify([model_meta.model_id for model_meta in storage.MODELS if model_meta.species == species])


def model_get(model_id):
    try:
        return jsonify(model_to_dict(storage.get(model_id).model))
    except KeyError:
        return f"Unknown model {model_id}", 404


def model_get_modified(model_id):
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    try:
        model_meta = storage.get(model_id)
    except KeyError:
        return f"Unknown model {model_id}", 404

    # Make a copy of the shared model instance for this request. It is not sufficient to use the cobra model context
    # manager here, as long as we're using async gunicorn workers and app state can be shared between requests.
    # This is an expensive operation, it can take a few seconds for large models.
    model = model_meta.model.copy()

    try:
        apply_operations(model, request.json['operations'])
    except KeyError:
        return "Missing field 'operations'", 400

    return jsonify(model_to_dict(model))


def model_modify(model_id):
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    try:
        model_meta = storage.get(model_id)
    except KeyError:
        return f"Unknown model '{model_id}'", 400

    # Make a copy of the shared model instance for this request. It is not sufficient to use the cobra model context
    # manager here, as long as we're using async gunicorn workers and app state can be shared between requests.
    # This is an expensive operation, it can take a few seconds for large models.
    model = model_meta.model.copy()

    try:
        conditions = request.json['conditions']
    except KeyError:
        return "Missing field 'conditions'", 400

    # Build list of operations to perform on the model
    operations = []
    errors = []
    if 'medium' in conditions:
        operations_medium, errors_medium = adapt_from_medium(model, conditions['medium'])
        operations.extend(operations_medium)
        errors.extend(errors_medium)

    if 'genotype' in conditions:
        operations_genotype, errors_genotype = adapt_from_genotype(model, conditions['genotype'])
        operations.extend(operations_genotype)
        errors.extend(errors_genotype)

    if 'measurements' in conditions:
        operations_measurements, errors_measurements = adapt_from_measurements(model, conditions['measurements'])
        operations.extend(operations_measurements)
        errors.extend(errors_measurements)

    if errors:
        # If any errors occured during modifications, discard generated operations and return the error messages to the
        # client for follow-up
        return jsonify({'errors': errors})
    else:
        delta_id = deltas.save(model.id, conditions, operations)
        return jsonify({'id': delta_id, 'operations': operations})


def model_simulate(model_id):
    if not request.is_json:
        return "Non-JSON request content is not supported", 415

    try:
        model_meta = storage.get(model_id)
    except KeyError:
        return f"Unknown model {model_id}", 404

    # Make a copy of the shared model instance for this request. It is not sufficient to use the cobra model context
    # manager here, as long as we're using async gunicorn workers and app state can be shared between requests.
    # This is an expensive operation, it can take a few seconds for large models.
    model = model_meta.model.copy()

    operations = []

    if 'delta_id' in request.json:
        delta_id = request.json['delta_id']
        try:
            operations.extend(deltas.load_from_key(delta_id))
        except KeyError:
            return f"Cannot find delta id '{delta_id}'", 404

    if 'operations' in request.json:
        operations.extend(request.json['operations'])

    apply_operations(model, operations)

    # Parse solver request args and set defaults
    method = request.json.get('method', 'fba')
    objective_id = request.json.get('objective')
    objective_direction = request.json.get('objective_direction')

    flux_distribution, growth_rate = simulate(model, method, objective_id, objective_direction)
    return jsonify({'flux_distribution': flux_distribution, 'growth_rate': growth_rate})


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
