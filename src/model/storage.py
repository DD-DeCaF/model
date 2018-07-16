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

from cobra.io import load_json_model, read_sbml_model
from redis import Redis

from model import constants
from model.app import app
from model.operations import restore_changes


logger = logging.getLogger(__name__)
redis = Redis(app.config['REDIS_ADDR'], app.config['REDIS_PORT'])


class ModelMeta:
    """A metabolic model, with metadata and an internal (lazy-loaded in development) cobrapy model instance."""

    def __init__(self, model_id, species, namespace, growth_rate_reaction, file_format='sbml', simulatable=True):
        self.model_id = model_id
        self.species = species
        self.namespace = namespace
        self.growth_rate_reaction = growth_rate_reaction
        self.file_format = file_format
        self.simulatable = simulatable

        # Preload the model into memory only in the following environments
        if app.config['ENVIRONMENT'] in ('production', 'staging'):
            self._load()

    @property
    def model(self):
        if not hasattr(self, '_model'):
            self._load()
        return self._model

    def _load(self):
        logger.info(f"Loading model {self.model_id} from SBML file")
        if self.file_format == 'sbml':
            self._model = read_sbml_model(f"data/models/{self.model_id}.sbml.gz")
        elif self.file_format == 'json':
            self._model = load_json_model(f"data/models/{self.model_id}.json")

        # Use cplex solver. Not necessary for the universal models.
        if self.simulatable:
            self._model.solver = 'cplex'
            self._model.notes['namespace'] = self.namespace


MODELS = [
    ModelMeta('iJO1366', 'ECOLX', 'bigg', 'BIOMASS_Ec_iJO1366_core_53p95M'),
    ModelMeta('iMM904', 'YEAST', 'bigg', 'BIOMASS_SC5_notrace'),
    ModelMeta('iMM1415', 'CRIGR', 'bigg', 'BIOMASS_mm_1_no_glygln'),
    ModelMeta('iNJ661', 'CORGT', 'bigg', 'BIOMASS_Mtb_9_60atp'),
    ModelMeta('iJN746', 'PSEPU', 'bigg', 'BIOMASS_KT_TEMP'),
    ModelMeta('e_coli_core', 'ECOLX', 'bigg', 'BIOMASS_Ecoli_core_w_GAM'),
    ModelMeta('ecYeast7', 'YEAST', 'yeast7', 'r_2111'),
    ModelMeta('ecYeast7_proteomics', 'YEAST', 'yeast7', 'r_2111'),
    ModelMeta('metanetx_universal_model_bigg', None, None, None, 'json', False),
]


def get(model_id):
    """Return the meta instance for the given model id"""
    for model_meta in MODELS:
        if model_meta.model_id == model_id:
            return model_meta
    else:
        raise KeyError(f"No model with id '{model_id}'")


def save_changes(model, message):
    """Save adapted changes in cache based on the given model and message"""
    changes_key = _changes_key(model.id, message)
    value = json.dumps({'model': model.id, 'changes': model.notes.get('changes', constants.get_empty_changes())})
    redis.set(changes_key, value)
    logger.info(f"Stored changes for {model.id} as '{changes_key}'")
    return changes_key


def restore_from_message(model_id, message):
    """Restore model with modifications from cache based on the given model and message"""
    return restore_from_key(_changes_key(model_id, message))


def restore_from_key(changes_key):
    """Restore model with modifications from cache based on the given unique cache key"""
    changes = redis.get(changes_key)
    if changes is None:
        logger.debug(f"Changes '{changes_key}' not found in cache")
        raise KeyError(f"Changes '{changes_key}' not found in cache")

    changes = json.loads(changes.decode('utf-8'))
    model = get(changes['model']).model.copy()
    model = restore_changes(model, changes['changes'])
    model.notes['changes'] = changes['changes']
    return model


def _changes_key(model_id, message):
    """Generate unique string used for cache key based on the model id and full message"""
    d = {k: message.get(k, []) for k in constants.MESSAGE_HASH_KEYS}
    d['model_id'] = model_id
    return hashlib.sha224(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()
