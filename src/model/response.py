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

import json
import logging
import os
import time

from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.flux_analysis.moma import add_moma
from cobra.io.dict import model_to_dict

from model import constants, storage
from model.operations import is_dummy, phase_plane_to_dict


logger = logging.getLogger(__name__)

METHODS = [
    'fba',
    'pfba',
    'fva',
    'pfba-fva',
    'moma',
    'lmoma',
]


def solve(model, method):
    if method not in METHODS:
        raise ValueError(f"Unsupported solver method '{method}'")

    try:
        if method == 'fba':
            solution = model.optimize()
        elif method == 'pfba':
            solution = pfba(model)
        elif method == 'fva':
            # FIXME: accept list of relevant fva reactions to calculate
            solution = flux_variability_analysis(model)
        elif method == 'pfba-fva':
            # FIXME: accept list of relevant fva reactions to calculate
            solution = flux_variability_analysis(model, fraction_of_optimum=1, pfba_factor=1.05)
        elif method in ('moma', 'lmoma'):
            reference = pfba(model)
            with model:
                add_moma(model, solution=reference, linear=(method == 'lmoma'))
                solution = model.optimize()
    except OptimizationError:
        flux = {}
        growth_rate = 0.0

        # For non-fva methods, do return the measured fluxes despite infesability
        if method not in ('fva', 'pfba-fva'):
            if 'measured' in changes:
                ids_measured_reactions = set(rxn['id'] for rxn in changes['measured']['reactions'])
                flux = {
                    rxn.id: (rxn.upper_bound + rxn.lower_bound) / 2
                    for rxn in model.reactions if rxn.id in ids_measured_reactions
                }
        return (flux, growth_rate)
    else:
        if method in ('fba', 'pfba', 'moma', 'lmoma'):
            flux = solution.fluxes.to_dict()
            growth_rate = flux[storage.get(model.id).growth_rate_reaction]
        elif method in ('fva', 'pfba-fva'):
            df = solution.rename(index=str, columns={'maximum': 'upper_bound', 'minimum': 'lower_bound'})
            for key in ['lower_bound', 'upper_bound']:
                df[key] = df[key].astype('float')
            flux = df.T.to_dict()
            growth_rate = flux[storage.get(model.id).growth_rate_reaction]['upper_bound']
        return (flux, growth_rate)


def respond(model, message, mutated_model_id=None):
    with model:
        if constants.OBJECTIVE in message:
            model.objective = model.reactions.get_by_id(message[constants.OBJECTIVE])
            if constants.OBJECTIVE_DIRECTION in message:
                model.objective.direction = message[constants.OBJECTIVE_DIRECTION]
        method = message.get(constants.SIMULATION_METHOD, 'fba')
        flux_distribution, growth_rate = solve(model, method)

    # Is it ok, to leave to-return optional?
    to_return = message.get('to-return', [
        constants.FLUXES,
        constants.TMY,
        constants.MODEL,
        constants.GROWTH_RATE,
        constants.REMOVED_REACTIONS,
        constants.MEASURED_REACTIONS,
        constants.ADDED_REACTIONS,
        constants.MISSING_MEASURED_REACTIONS,
    ])

    # Add requested results to the response
    response = {}

    if constants.FLUXES in to_return:
        response[constants.FLUXES] = flux_distribution

    if constants.TMY in to_return:
        objectives = message.get(constants.TMY_OBJECTIVES, [])
        response[constants.TMY] = {key: phase_plane_to_dict(model, key) for key in objectives}

    if constants.MODEL in to_return:
        response[constants.MODEL] = model_to_dict(model)

    if constants.GROWTH_RATE in to_return:
        response[constants.GROWTH_RATE] = growth_rate

    if constants.REMOVED_REACTIONS in to_return:
        changes = model.notes.get('changes', constants.get_empty_changes())
        reactions = list(set([i['id'] for i in changes['removed']['reactions'] ]))
        response[constants.REMOVED_REACTIONS] = reactions

    if constants.MEASURED_REACTIONS in to_return:
        changes = model.notes.get('changes', constants.get_empty_changes())
        reactions = list(set([i['id'] for i in changes['measured']['reactions']] ))
        response[constants.MEASURED_REACTIONS] = reactions

    if constants.ADDED_REACTIONS in to_return:
        changes = model.notes.get('changes', constants.get_empty_changes())
        reactions = list(set([i['id'] for i in changes['added']['reactions'] if not is_dummy(i['id'])]))
        response[constants.ADDED_REACTIONS] = reactions

    if constants.MISSING_MEASURED_REACTIONS in to_return:
        changes = model.notes.get('changes', constants.get_empty_changes())
        reactions = list(set([i['id'] for i in changes['measured-missing']['reactions']]))
        response[constants.MISSING_MEASURED_REACTIONS] = reactions

    if mutated_model_id:
        response['model-id'] = mutated_model_id

    return response
