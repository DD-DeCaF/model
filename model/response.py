import jsonpatch
import time

from cobra.exceptions import OptimizationError
from cobra.flux_analysis import pfba
from cobra.io.dict import model_to_dict

import model.consts as consts
from model.logger import logger
import model.utils as utils
from model.model_operations import is_dummy, phase_plane_to_dict

class Response(object):
    def __init__(self, model, message, wild_type_model=None):
        self.model = model
        self.message = message
        self.method_name = message.get(consts.SIMULATION_METHOD, 'fba')
        self.wild_type_model = wild_type_model
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
                self.growth = self.flux[consts.MODEL_GROWTH_RATE[model.id]]['upper_bound']
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
                self.growth = self.flux[consts.MODEL_GROWTH_RATE[model.id]]

    def solve_fva(self):
        fva_reactions = None
        if consts.MAP in self.message:
            reaction_ids = utils.map_reactions_list('{0}/{1}/{1}.{2}.json'.format(
                consts.MAPS_DIR,
                self.model.id,
                self.message[consts.MAP]
            ))
            if reaction_ids:
                reactions = [i for i in reaction_ids
                             if self.model.reactions.has_id(i)]
                fva_reactions = list(set(
                    reactions + [consts.MODEL_GROWTH_RATE[self.model.id]]
                ))
        return consts.METHODS[self.method_name](
            self.model,
            reactions=fva_reactions
        )

    def solve(self):
        t = time.time()
        if self.method_name in {'lmoma', 'moma'}:
            solution = consts.METHODS[self.method_name](self.model, reference=pfba(self.model))
        else:
            solution = consts.METHODS[self.method_name](self.model)
        logger.info('Model solved with method %s in %s sec', self.method_name, time.time() - t)
        return solution

    def model_json(self):
        if self.wild_type_model is not None:
            return jsonpatch.make_patch(model_to_dict(self.wild_type_model), model_to_dict(self.model)).patch
        else:
            return model_to_dict(self.model)

    def fluxes(self):
        return self.flux

    def theoretical_maximum_yield(self):
        objectives = self.message.get(consts.OBJECTIVES, [])
        res = {key: phase_plane_to_dict(self.model, key) for key in objectives}
        logger.info(res)
        return res

    def growth_rate(self):
        return self.growth

    def measured_missing_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', consts.get_empty_changes())['measured-missing']['reactions']]
        ))

    def removed_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', consts.get_empty_changes())['removed']['reactions']]
        ))

    def added_reactions(self):
        return list(set(
            [i['id'] for i in self.model.notes.get('changes', consts.get_empty_changes())['added']['reactions']
             if not is_dummy(i['id'])]
        ))

    def full_response(self):
        ret = {}
        for key, func in consts.RETURN_FUNCTIONS.items():
            ret[key] = getattr(self, consts.RETURN_FUNCTIONS[key])()
        return ret


async def respond(model, message=None, mutated_model_id=None, wild_type_model=None):
    message = message if message is not None else {}
    t = time.time()
    response = Response(model, message, wild_type_model)
    # Is it ok, to leave to-return optional?
    to_return = message.get('to-return', None)
    if to_return is not None:
        result = {}
        for key in message['to-return']:
            result[key] = getattr(response, consts.RETURN_FUNCTIONS[key])()
    else:
        result = response.full_response()
    if mutated_model_id:
        result['model-id'] = mutated_model_id
    if consts.REQUEST_ID in message:
        result[consts.REQUEST_ID] = message[consts.REQUEST_ID]
    logger.info('Response for %s is ready in %s sec', message, time.time() - t)
    return result
