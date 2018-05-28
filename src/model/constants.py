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
import os
import re
import time
from itertools import chain

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
REQUEST_ID = 'request-id'
REMOVED_REACTIONS = 'removed-reactions'
ADDED_REACTIONS = 'added-reactions'
MISSING_MEASURED_REACTIONS = 'missing-measured-reactions'
MEASURED_REACTIONS = 'measured-reactions'
MESSAGE_HASH_KEYS = {GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS, REACTIONS_ADD, REACTIONS_KNOCKOUT, OBJECTIVE,
                     MEASURED_REACTIONS}


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

SPECIES_TO_MODEL = {
    'ECOLX': ['iJO1366', 'e_coli_core'],
    'YEAST': ['iMM904', 'ecYeast7', 'ecYeast7_proteomics'],
    'CRIGR': ['iMM1415'],
    'CORGT': ['iNJ661'],
    'PSEPU': ['iJN746'],
}

MODELS = frozenset(chain.from_iterable(SPECIES_TO_MODEL.values()))

MODEL_NAMESPACE = {
    'iJO1366': 'bigg',
    'iMM904': 'bigg',
    'iMM1415': 'bigg',
    'iNJ661': 'bigg',
    'iJN746': 'bigg',
    'e_coli_core': 'bigg',
    'ecYeast7': 'yeast7',
    'ecYeast7_proteomics': 'yeast7',
}

MODEL_GROWTH_RATE = {
    'iJO1366': 'BIOMASS_Ec_iJO1366_core_53p95M',
    'iMM904': 'BIOMASS_SC5_notrace',
    'iMM1415': 'BIOMASS_mm_1_no_glygln',
    'iNJ661': 'BIOMASS_Mtb_9_60atp',
    'iJN746': 'BIOMASS_KT_TEMP',
    'e_coli_core': 'BIOMASS_Ecoli_core_w_GAM',
    'ecYeast7': 'r_2111',
    'ecYeast7_proteomics': 'r_2111',
}

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


def generate_map_dictionary():
    """Generate model-maps lookup depending on the folder structure

    :return: dict
    """
    result = {}
    for path, _, files in os.walk(MAPS_DIR):
        if files:
            result[path.replace(MAPS_DIR + '/', '')] = \
                sorted([re.match(r".*\.(.+)\..*", f).group(1) for f in files])
    return result


MAP_DICTIONARY = generate_map_dictionary()
