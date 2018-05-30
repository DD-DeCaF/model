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

import jsonpatch
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import pfba
from cobra.io.dict import model_to_dict

import model.constants as constants
from model.operations import is_dummy, phase_plane_to_dict
from model.storage import Models


logger = logging.getLogger(__name__)


def map_reactions_list(map_path):
    """Extract reaction ids from map for FVA optimization

    :param map_path: string
    :return: list of strings
    """
    if not os.path.isfile(map_path):
        return []
    with open(map_path) as fp:
        return [i['bigg_id'] for i in json.load(fp)[1]['reactions'].values()]


class Response(object):
    def __init__(self, model, message, wild_type_model_id=None):
        self.model = model
        with self.model:
            if constants.OBJECTIVE in message:
                self.model.objective = self.model.reactions.get_by_id(message[constants.OBJECTIVE])
            self.message = message
            self.method_name = message.get(constants.SIMULATION_METHOD, 'fba')
            self.wild_type_model_id = wild_type_model_id
            if self.method_name in {'fva', 'pfba-fva'}:
                try:
                    solution = self.solve_fva()
                except OptimizationError:
                    logger.info('infeasible model for fva')
                    self.flux = {}
                    self.growth = 0.0
                else:
                    df = solution.rename(index=str, columns={'maximum': 'upper_bound', 'minimum': 'lower_bound'})
                    for key in ['lower_bound', 'upper_bound']:
                        df[key] = df[key].astype('float')
                    self.flux = df.T.to_dict()
                    self.growth = self.flux[constants.MODEL_GROWTH_RATE[model.id]]['upper_bound']
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
                    self.growth = self.flux[constants.MODEL_GROWTH_RATE[model.id]]

    def solve_fva(self):
        fva_reactions = None
        if constants.MAP in self.message:
            reaction_ids = map_reactions_list('{0}/{1}/{1}.{2}.json'.format(
                constants.MAPS_DIR,
                self.model.id,
                self.message[constants.MAP]
            ))
            if reaction_ids:
                reactions = [i for i in reaction_ids
                             if self.model.reactions.has_id(i)]
                fva_reactions = list(set(
                    reactions + [constants.MODEL_GROWTH_RATE[self.model.id]]
                ))
        return constants.METHODS[self.method_name](
            self.model,
            reactions=fva_reactions
        )

    def solve(self):
        t = time.time()
        if self.method_name in {'lmoma', 'moma'}:
            solution = constants.METHODS[self.method_name](self.model, reference=pfba(self.model))
        else:
            solution = constants.METHODS[self.method_name](self.model)
        logger.info('Model solved with method %s in %s sec', self.method_name, time.time() - t)
        return solution

    def model_json(self):
        if self.wild_type_model_id is not None:
            return jsonpatch.make_patch(Models.get_dict(self.wild_type_model_id), model_to_dict(self.model)).patch
        else:
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


async def respond(model, message=None, mutated_model_id=None, wild_type_model_id=None):
    message = message if message is not None else {}
    t = time.time()
    response = Response(model, message, wild_type_model_id)
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
    if constants.REQUEST_ID in message:
        result[constants.REQUEST_ID] = message[constants.REQUEST_ID]
    logger.info('Response for %s is ready in %s sec', message, time.time() - t)
    return result
