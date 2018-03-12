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
from copy import deepcopy

import aiohttp
import jsonpatch
from aiohttp.test_utils import AioHTTPTestCase, unittest_run_loop

from model.app import get_app


logging.disable(logging.CRITICAL)

MEASUREMENTS = [{'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0],
                 'type': 'compound'},
                {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9],
                 'type': 'compound'},
                {'id': 'PFK', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'},
                {'id': 'BAD_ID', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'}]
MESSAGE_FLUXES = {'to-return': ['fluxes'], 'measurements': MEASUREMENTS}
MESSAGE_FLUXES_INFEASIBLE = {'to-return': ['fluxes'], 'measurements': [
    {'id': 'ATPM', 'measurements': [100, 100], 'type': 'reaction', 'db_name': 'bigg.reaction'}]}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy', 'model'], 'theoretical-objectives': ['chebi:17790'], 'request-id': 'requestid'}
MESSAGE_MODIFY = {
    'simulation-method': 'pfba',
    'reactions-add': [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR83321', 'metabolites': None},
        {'id': 'MagicCarrot', 'metabolites': {'glc__D_c': -1, 'caro_c': 1}}
    ],
    'to-return': ['tmy', 'fluxes', 'model', 'added-reactions', 'removed-reactions'],
    'theoretical-objectives': ['chebi:17790', 'chebi:17579'],
    'genotype-changes': ['+Aac', '-pta'],
    'measurements': MEASUREMENTS,
}

MESSAGE_DIFFERENT_OBJECTIVE = {'to-return': ['fluxes'], 'objective': 'EX_etoh_e'}

MODELS_URL = '/models/{}'
V1_MODELS_URL = '/v1/models/{}'
MODEL_OPTIONS_URL = '/model-options/{}'
MAPS_URL = '/maps'
WS_URL = '/wsmodels/{}'


class EndToEndTestCase(AioHTTPTestCase):
    async def get_application(self):
        return get_app()

    @unittest_run_loop
    async def test_map_success(self):
        params = {'model': 'e_coli_core', 'map': 'Core metabolism'}
        response = await self.client.get('/map', params=params)
        response.raise_for_status()

    @unittest_run_loop
    async def test_map_traversal_attempt(self):
        response = await self.client.get('/map', params={'model': '../../e_coli_core', 'map': 'Core metabolism'})
        assert response.status == 400

    @unittest_run_loop
    async def test_map_not_found(self):
        response = await self.client.get('/map', params={'model': 'e_coli_core', 'map': 'Non existing'})
        assert response.status == 404

    @unittest_run_loop
    async def test_model_info(self):
        response = await self.client.get('/v1/model-info/{}'.format('e_coli_core'))
        assert (await response.json())['medium'] == [
            {'id': 'EX_co2_e', 'name': 'CO2'},
            {'id': 'EX_glc__D_e', 'name': 'D-Glucose'},
            {'id': 'EX_h2o_e', 'name': 'H2O'},
            {'id': 'EX_h_e', 'name': 'H+'},
            {'id': 'EX_nh4_e', 'name': 'Ammonia'},
            {'id': 'EX_o2_e', 'name': 'O2'},
            {'id': 'EX_pi_e', 'name': 'Phosphate'}
        ]

    @unittest_run_loop
    async def test_http(self):
        response = await self.client.get(MODEL_OPTIONS_URL.format('ECOLX'))
        response.raise_for_status()
        response = await self.client.get(MAPS_URL)
        response.raise_for_status()
        for query, message in {'modify': MESSAGE_MODIFY, 'fluxes': MESSAGE_FLUXES}.items():
            response = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': message})
            response.raise_for_status()
            model = await response.json()
            assert set(model.keys()) == set(message['to-return'] + ['model-id'])
            assert model['fluxes']['EX_glc__D_e'] == -9.0
            assert abs(model['fluxes']['EX_etoh_e'] - 4.64) < 0.001  # lower bound
            assert model['fluxes']['PFK'] == 5
            if query == 'modify':
                tmy = model['tmy']
                changes = model['model']['notes']['changes']
                assert sum(tmy['chebi:17579']['DM_caro_e']) > 10.
                assert 'PTA2' in {rxn['id'] for rxn in changes['removed']['reactions']}
                assert 'EX_etoh_e' in {rxn['id'] for rxn in changes['measured']['reactions']}
                assert 'PFK' in {rxn['id'] for rxn in changes['measured']['reactions']}
                assert 'b2297' in {rxn['id'] for rxn in changes['removed']['genes']}
                assert 'BAD_ID' in {rxn['id'] for rxn in changes['measured-missing']['reactions']}
        response = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': MESSAGE_FLUXES_INFEASIBLE})
        model = await response.json()
        assert model['fluxes']['ATPM'] == 100
        response = await self.client.post(MODELS_URL.format('wrong_id'), json={'message': {}})
        assert response.status == 404
        response = await self.client.post(MODELS_URL.format('iJO1366'), json={})
        assert response.status == 400
        response_etoh = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': MESSAGE_DIFFERENT_OBJECTIVE})
        response_etoh.raise_for_status()
        model = await response_etoh.json()
        assert abs(model['fluxes']['EX_etoh_e']) - 20.0 < 0.001
        MESSAGE_DIFFERENT_OBJECTIVE.pop('objective')
        response_etoh = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': MESSAGE_DIFFERENT_OBJECTIVE})
        response_etoh.raise_for_status()
        model = await response_etoh.json()
        assert abs(model['fluxes']['EX_etoh_e']) < 0.001

    @unittest_run_loop
    async def test_websocket(self):
        response = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': MESSAGE_MODIFY})
        response.raise_for_status()
        model_id = (await response.json())['model-id']
        ws = await self.client.ws_connect(WS_URL.format(model_id))
        ws.send_json(MESSAGE_TMY_FLUXES)
        async for msg in ws:
            if msg.type == aiohttp.WSMsgType.TEXT:
                msg_content = msg.json()
                assert msg_content['fluxes']
                assert msg_content['tmy']
                assert isinstance(msg_content['model'], dict)
            elif msg.type == aiohttp.WSMsgType.ERROR:
                raise ws.exception()
            await ws.close()
            break

    @unittest_run_loop
    async def test_websocket_v1(self):
        response = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': MESSAGE_MODIFY})
        response.raise_for_status()
        model_id = (await response.json())['model-id']
        ws = await self.client.ws_connect('/v1{}'.format(WS_URL.format(model_id)))
        ws.send_json(MESSAGE_TMY_FLUXES)
        async for msg in ws:
            if msg.type == aiohttp.WSMsgType.TEXT:
                msg_content = msg.json()
                assert msg_content['fluxes']
                assert msg_content['tmy']
                assert isinstance(msg_content['model'], list)
            elif msg.type == aiohttp.WSMsgType.ERROR:
                raise ws.exception()
            await ws.close()
            break

    @unittest_run_loop
    async def test_diff_model_api(self):
        original_message = deepcopy(MESSAGE_MODIFY)
        original_message['to-return'] = ["fluxes", "model", "added-reactions", "removed-reactions"]
        original_response = await self.client.post(MODELS_URL.format('iJO1366'), json={'message': original_message})
        try:
            original_response.raise_for_status()
        except Exception as ex:
            print('Original resp: ', await original_response.text())
            raise ex
        original_model = (await original_response.json())['model']

        wild_type_model_response = await self.client.get(V1_MODELS_URL.format('iJO1366'))
        wild_type_model_response.raise_for_status()
        wild_type_model = await wild_type_model_response.json()

        diff_response = await self.client.post(V1_MODELS_URL.format('iJO1366'), json={'message': original_message})
        diff_response.raise_for_status()
        patch = jsonpatch.JsonPatch((await diff_response.json())['model'])
        result = patch.apply(wild_type_model)

        del result['notes']['changes']
        del original_model['notes']['changes']
        assert result == original_model
