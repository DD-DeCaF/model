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
from cobra.flux_analysis import pfba
from cobra.io.dict import model_to_dict

from model import constants, storage
from model.operations import is_dummy, phase_plane_to_dict


logger = logging.getLogger(__name__)


class Response(object):
    def __init__(self, model, message):
        self.model = model
        with self.model:
            if constants.OBJECTIVE in message:
                self.model.objective = self.model.reactions.get_by_id(message[constants.OBJECTIVE])
                if constants.OBJECTIVE_DIRECTION in message:
                    self.model.objective.direction = message[constants.OBJECTIVE_DIRECTION]
            self.message = message
            self.method_name = message.get(constants.SIMULATION_METHOD, 'fba')
            if self.method_name in {'fva', 'pfba-fva'}:
                try:
                    # FIXME: accept list of reactions to actually solve for fva
                    solution = constants.METHODS[self.method_name](self.model)
                except OptimizationError:
                    logger.info('infeasible model for fva')
                    self.flux = {}
                    self.growth = 0.0
                else:
                    df = solution.rename(index=str, columns={'maximum': 'upper_bound', 'minimum': 'lower_bound'})
                    for key in ['lower_bound', 'upper_bound']:
                        df[key] = df[key].astype('float')
                    self.flux = df.T.to_dict()
                    self.growth = self.flux[storage.get(model.id).growth_rate_reaction]['upper_bound']
            else:
                try:
                    solution = self.solve()
                except OptimizationError:
                    logger.info('infeasible model, returning measured fluxes only')
                    changes = model.notes['changes']
                    self.flux = {}
                    self.growth = 0.0
                    if 'measured' in changes:
                        ids_measured_reactions = set(rxn['id'] for rxn in changes['measured']['reactions'])
                        self.flux = dict((rxn.id, (rxn.upper_bound + rxn.lower_bound) / 2)
                                         for rxn in model.reactions if rxn.id in ids_measured_reactions)
                else:
                    self.flux = solution.fluxes.to_dict()
                    self.growth = self.flux[storage.get(model.id).growth_rate_reaction]

    def solve(self):
        t = time.time()
        if self.method_name in {'lmoma', 'moma'}:
            solution = constants.METHODS[self.method_name](self.model, reference=pfba(self.model))
        else:
            solution = constants.METHODS[self.method_name](self.model)
        logger.info('Model solved with method %s in %s sec', self.method_name, time.time() - t)
        return solution

    def model_json(self):
        return model_to_dict(self.model)

    def fluxes(self):
        return self.flux

    def theoretical_maximum_yield(self):
        objectives = self.message.get(constants.TMY_OBJECTIVES, [])
        res = {key: phase_plane_to_dict(self.model, key) for key in objectives}
        logger.info(res)
        return res

    def growth_rate(self):
        return self.growth

    def measured_missing_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes',
                                                   constants.get_empty_changes())['measured-missing']['reactions']]
        ))

    def removed_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', constants.get_empty_changes())['removed']['reactions']]
        ))

    def added_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', constants.get_empty_changes())['added']['reactions']
             if not is_dummy(i['id'])]
        ))

    def measured_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', constants.get_empty_changes())['measured']['reactions']]
        ))

    def full_response(self):
        ret = {}
        for key, func in constants.RETURN_FUNCTIONS.items():
            ret[key] = getattr(self, constants.RETURN_FUNCTIONS[key])()
        return ret


def respond(model, message=None, mutated_model_id=None):
    message = message if message is not None else {}
    t = time.time()
    response = Response(model, message)
    # Is it ok, to leave to-return optional?
    to_return = message.get('to-return', None)
    if to_return is not None:
        result = {}
        for key in message['to-return']:
            result[key] = getattr(response, constants.RETURN_FUNCTIONS[key])()
    else:
        result = response.full_response()
    if mutated_model_id:
        result['model-id'] = mutated_model_id
    logger.info('Response for %s is ready in %s sec', message, time.time() - t)
    return result
