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

from cobra.io.dict import model_from_dict
import requests

from model.app import app


logger = logging.getLogger(__name__)


class ModelWrapper:
    """A wrapper for a cobrapy model with some additional metadata."""

    def __init__(self, model, organism_id, biomass_reaction):
        self.model = model
        # Use the cplex solver for performance
        self.model.solver = 'cplex'
        self.organism_id = organism_id
        self.biomass_reaction = biomass_reaction


# Keep all loaded models in memory in this dictionary, keyed by our internal
# model storage primary key id.
_MODELS = {}


def get(model_id):
    """Return a ModelWrapper instance for the given model id"""
    if model_id not in _MODELS:
        _load_model(model_id)
    return _MODELS[model_id]


def _load_model(model_id):
    logger.debug(f"Requesting model {model_id} from the model warehouse")
    response = requests.get(f"{app.config['MODEL_WAREHOUSE_API']}/models/{model_id}")

    if response.status_code == 404:
        raise KeyError(f"No model with id {model_id}")
    response.raise_for_status()

    logger.debug(f"Deserializing received model with cobrapy")
    model_data = response.json()
    _MODELS[model_id] = ModelWrapper(
        model_from_dict(model_data['model_serialized']),
        model_data['organism_id'],
        model_data['default_biomass_reaction'],
    )


# Preload all models in production/staging environments
if app.config['ENVIRONMENT'] in ('production', 'staging'):
    logger.info(f"Preloading all public models (this may take some time)")
    response = requests.get(f"{app.config['MODEL_WAREHOUSE_API']}/models")
    response.raise_for_status()
    for model in response.json():
        _load_model(model['id'])
    logger.info(f"Done preloading {len(response.json())} models")
