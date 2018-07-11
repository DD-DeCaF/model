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
import time

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.flux_analysis.moma import add_moma


logger = logging.getLogger(__name__)


def pfba_fva(model, reactions=None):
    return flux_variability_analysis(
        model,
        fraction_of_optimum=1,
        pfba_factor=1.05,
        reactions_list=reactions
    )


def moma(model, reference, linear=False):
    with model:
        start_time = time.time()
        add_moma(model, solution=reference, linear=linear)
        logger.info('moma addition finished in %s s', time.time() - start_time)
        start_time = time.time()
        solution = model.optimize()
        logger.info('moma optimization finished in %s s', time.time() - start_time)
    return solution


METHODS = {
    'fba': lambda model: model.optimize(),
    'pfba': pfba,
    'fva': flux_variability_analysis,
    'pfba-fva': pfba_fva,
    'moma': moma,
    'lmoma': lambda model, reference: moma(model, reference, linear=True),
}

GENOTYPE_CHANGES = 'genotype-changes'
MEDIUM = 'medium'
MEASUREMENTS = 'measurements'
SIMULATION_METHOD = 'simulation-method'
MAP = 'map'
REACTIONS_KNOCKOUT = 'reactions-knockout'
REACTIONS_ADD = 'reactions-add'
MODEL = 'model'
FLUXES = 'fluxes'
GROWTH_RATE = 'growth-rate'
TMY = 'tmy'
TMY_OBJECTIVES = 'theoretical-objectives'
OBJECTIVE = 'objective'
REMOVED_REACTIONS = 'removed-reactions'
ADDED_REACTIONS = 'added-reactions'
MISSING_MEASURED_REACTIONS = 'missing-measured-reactions'
MEASURED_REACTIONS = 'measured-reactions'
OBJECTIVE_DIRECTION = 'objective-direction'

MESSAGE_HASH_KEYS = {GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS, REACTIONS_ADD, REACTIONS_KNOCKOUT, OBJECTIVE,
                     OBJECTIVE_DIRECTION, MEASURED_REACTIONS}


def get_empty_changes():
    return {
        'added': {
            'reactions': [],
            'metabolites': [],
        },
        'removed': {
            'genes': [],
            'reactions': [],
        },
        'measured': {
            'genes': [],
            'reactions': []
        },
        'measured-missing': {
            'genes': [],
            'reactions': []
        },
        'measured-medium': {
            'reactions': []
        }
    }


MAPS_DIR = 'data/maps'

REQUEST_KEYS = [GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS]

RETURN_FUNCTIONS = {
    FLUXES: 'fluxes',
    TMY: 'theoretical_maximum_yield',
    MODEL: 'model_json',
    GROWTH_RATE: 'growth_rate',
    REMOVED_REACTIONS: 'removed_reactions',
    MEASURED_REACTIONS: 'measured_reactions',
    ADDED_REACTIONS: 'added_reactions',
    MISSING_MEASURED_REACTIONS: 'measured_missing_reactions',
}
