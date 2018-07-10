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

from copy import deepcopy

import jsonpatch
import pytest


MEASUREMENTS = [{'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0],
                 'type': 'compound'},
                {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9],
                 'type': 'compound'},
                {'id': 'PFK', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'},
                {'id': 'BAD_ID', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'}]
MESSAGE_FLUXES = {'to-return': ['fluxes'], 'measurements': MEASUREMENTS}
MESSAGE_FLUXES_INFEASIBLE = {'to-return': ['fluxes'], 'measurements': [
    {'id': 'ATPM', 'measurements': [100, 100], 'type': 'reaction', 'db_name': 'bigg.reaction'}]}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy', 'model'], 'theoretical-objectives': ['chebi:17790'],
                      'request-id': 'requestid'}
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


def test_model_info(client):
    response = client.get('/v1/model-info/{}'.format('e_coli_core'))
    assert response.json['medium'] == [
        {'id': 'EX_co2_e', 'name': 'CO2'},
        {'id': 'EX_glc__D_e', 'name': 'D-Glucose'},
        {'id': 'EX_h2o_e', 'name': 'H2O'},
        {'id': 'EX_h_e', 'name': 'H+'},
        {'id': 'EX_nh4_e', 'name': 'Ammonia'},
        {'id': 'EX_o2_e', 'name': 'O2'},
        {'id': 'EX_pi_e', 'name': 'Phosphate'}
    ]


def test_http(client):
    response = client.get(MODEL_OPTIONS_URL.format('ECOLX'))
    assert response.status_code == 200
    for query, message in {'modify': MESSAGE_MODIFY, 'fluxes': MESSAGE_FLUXES}.items():
        response = client.post(MODELS_URL.format('iJO1366'), json={'message': message})
        assert response.status_code == 200
        model = response.json
        assert set(model.keys()) == set(message['to-return'] + ['model-id'])
        assert model['fluxes']['EX_glc__D_e'] == -9.0
        assert abs(model['fluxes']['EX_etoh_e'] - 4.64) < 0.001  # lower bound
        assert model['fluxes']['PFK'] == pytest.approx(5)
        if query == 'modify':
            tmy = model['tmy']
            changes = model['model']['notes']['changes']
            assert sum(tmy['chebi:17579']['DM_caro_e']) > 10.
            assert 'PTA2' in {rxn['id'] for rxn in changes['removed']['reactions']}
            assert 'EX_etoh_e' in {rxn['id'] for rxn in changes['measured']['reactions']}
            assert 'PFK' in {rxn['id'] for rxn in changes['measured']['reactions']}
            assert 'b2297' in {rxn['id'] for rxn in changes['removed']['genes']}
            assert 'BAD_ID' in {rxn['id'] for rxn in changes['measured-missing']['reactions']}
    response = client.post(MODELS_URL.format('iJO1366'), json={'message': MESSAGE_FLUXES_INFEASIBLE})
    model = response.json
    assert model['fluxes']['ATPM'] == pytest.approx(100)
    response = client.post(MODELS_URL.format('wrong_id'), json={'message': {}})
    assert response.status_code == 404
    response = client.post(MODELS_URL.format('iJO1366'), json={})
    assert response.status_code == 400
    response_etoh = client.post(MODELS_URL.format('iJO1366'),
                                json={'message': MESSAGE_DIFFERENT_OBJECTIVE})
    assert response_etoh.status_code == 200
    model = response_etoh.json
    assert abs(model['fluxes']['EX_etoh_e']) - 20.0 < 0.001
    MESSAGE_DIFFERENT_OBJECTIVE.pop('objective')
    response_etoh = client.post(MODELS_URL.format('iJO1366'),
                                json={'message': MESSAGE_DIFFERENT_OBJECTIVE})
    assert response_etoh.status_code == 200
    model = response_etoh.json
    assert abs(model['fluxes']['EX_etoh_e']) < 0.001
