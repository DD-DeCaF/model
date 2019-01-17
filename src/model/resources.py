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

from cobra.io.dict import model_from_dict, model_to_dict
from flask import Response, abort, jsonify, request
from flask_apispec import use_kwargs
from flask_apispec.extension import FlaskApiSpec
from prometheus_client import CONTENT_TYPE_LATEST, CollectorRegistry, generate_latest
from prometheus_client.multiprocess import MultiProcessCollector

from model import storage
from model.exceptions import Forbidden, ModelNotFound, Unauthorized
from model.modeling.adapter import apply_genotype, apply_measurements, apply_medium
from model.modeling.operations import apply_operations
from model.modeling.simulations import simulate
from model.schemas import ModificationRequest, Operation, SimulationRequest


logger = logging.getLogger(__name__)


def init_app(app):
    """Register API resources on the provided Flask application."""
    app.add_url_rule('/healthz', view_func=healthz)
    app.add_url_rule('/metrics', view_func=metrics)

    app.add_url_rule('/models/<model_id>', view_func=model_get_modified, methods=['POST'])
    app.add_url_rule('/models/<model_id>/modify', view_func=model_modify, methods=['POST'])
    app.add_url_rule('/simulate', view_func=model_simulate, methods=['POST'])

    docs = FlaskApiSpec(app)
    docs.register(model_get_modified, endpoint=model_get_modified.__name__)
    docs.register(model_modify, endpoint=model_modify.__name__)
    docs.register(model_simulate, endpoint=model_simulate.__name__)


@use_kwargs(Operation(many=True))
def model_get_modified(model_id, operations):
    try:
        model_wrapper = storage.get(model_id)
    except Unauthorized as error:
        abort(401, error.message)
    except Forbidden as error:
        abort(403, error.message)
    except ModelNotFound as error:
        abort(404, error.message)

    # Use the context manager to undo all modifications to the shared model instance on completion.
    with model_wrapper.model as model:
        apply_operations(model, operations)
        return jsonify(model_to_dict(model))


@use_kwargs(ModificationRequest)
def model_modify(model_id, medium, genotype, measurements):
    if not request.is_json:
        abort(415, "Non-JSON request content is not supported")

    try:
        model_wrapper = storage.get(model_id)
    except Unauthorized as error:
        abort(401, error.message)
    except Forbidden as error:
        abort(403, error.message)
    except ModelNotFound as error:
        abort(404, error.message)

    # Use the context manager to undo all modifications to the shared model instance on completion.
    with model_wrapper.model as model:
        # Build list of operations to perform on the model
        operations = []
        errors = []
        if medium:
            operations_medium, errors_medium = apply_medium(model, medium)
            operations.extend(operations_medium)
            errors.extend(errors_medium)

        if genotype:
            operations_genotype, errors_genotype = apply_genotype(model, genotype)
            operations.extend(operations_genotype)
            errors.extend(errors_genotype)

        if measurements:
            operations_measurements, errors_measurements = apply_measurements(
                model,
                model_wrapper.biomass_reaction,
                measurements,
            )
            operations.extend(operations_measurements)
            errors.extend(errors_measurements)

        if errors:
            # If any errors occured during modifications, discard generated operations and return the error messages to
            # the client for follow-up
            return jsonify({'errors': errors}), 400
        else:
            return jsonify({'operations': operations})


@use_kwargs(SimulationRequest)
def model_simulate(model_id, model, biomass_reaction, method, objective_id, objective_direction, operations):
    if model_id:
        # Client provided an identifier to some model in the model storage service, retrieve it
        logger.debug(f"Simulating model by id {request.json['model_id']}")

        try:
            model_wrapper = storage.get(request.json['model_id'])
            biomass_reaction = model_wrapper.biomass_reaction
        except Unauthorized as error:
            abort(401, error.message)
        except Forbidden as error:
            abort(403, error.message)
        except ModelNotFound as error:
            abort(404, error.message)

        model = model_wrapper.model
    elif model:
        # Client provided a full custom model and biomass reaction
        logger.debug("Simulating ad-hoc model provided by client")

        try:
            model_dict = request.json['model']
            model = model_from_dict(model_dict)
        except Exception as error:
            logger.warning(f"Cobrapy could not deserialize provided model: {str(error)}", exc_info=True)
            logger.debug(f"Full serialized model: {model_dict}")
            abort(400, f"The provided model is not deserializable by cobrapy")
        if not biomass_reaction:
            abort(400, f"Missing field 'biomass_reaction'")
        if not model.reactions.has_id(biomass_reaction):
            abort(400, f"There is no biomass reaction with id '{biomass_reaction}' in the model")
    else:
        abort(400, f"Missing field 'model' or 'model_id'")

    # Use the context manager to undo all modifications to the shared model instance on completion.
    with model:
        apply_operations(model, operations)
        flux_distribution, growth_rate = simulate(model, biomass_reaction, method, objective_id, objective_direction)
        return jsonify({'flux_distribution': flux_distribution, 'growth_rate': growth_rate})


def metrics():
    return Response(generate_latest(MultiProcessCollector(CollectorRegistry())), mimetype=CONTENT_TYPE_LATEST)


def healthz():
    """
    HTTP endpoint for readiness probes.

    Return an empty response. This response will not be ready until the application has finished initializing, e.g.,
    preloading models, which takes a few minutes.
    """
    return ""
