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


logger = logging.getLogger(__name__)

GENOTYPE_CHANGES = 'genotype-changes'
MEDIUM = 'medium'
MEASUREMENTS = 'measurements'
SIMULATION_METHOD = 'simulation-method'
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


REQUEST_KEYS = [GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS]

