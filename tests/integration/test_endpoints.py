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

"""Test local HTTP endpoints"""

import pytest

from model import adapter
from model.ice_client import ICE


MEASUREMENTS = [{'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0],
                 'type': 'compound'},
                {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9],
                 'type': 'compound'},
                {'id': 'PFK', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'},
                {'id': 'BAD_ID', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'}]
MESSAGE_FLUXES = {'to-return': ['fluxes'], 'measurements': MEASUREMENTS}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy', 'model'], 'theoretical-objectives': ['chebi:17790']}
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


def test_model_info(client):
    response = client.get('/models/{}/medium'.format('e_coli_core'))
    assert response.status_code == 200
    assert response.json['medium'] == [
        {'id': 'EX_co2_e', 'name': 'CO2'},
        {'id': 'EX_glc__D_e', 'name': 'D-Glucose'},
        {'id': 'EX_h2o_e', 'name': 'H2O'},
        {'id': 'EX_h_e', 'name': 'H+'},
        {'id': 'EX_nh4_e', 'name': 'Ammonia'},
        {'id': 'EX_o2_e', 'name': 'O2'},
        {'id': 'EX_pi_e', 'name': 'Phosphate'}
    ]


def test_species(client):
    response = client.get("/species/ECOLX")
    assert response.status_code == 200
    assert len(response.json) > 0


def test_simulate_wrong_id(client):
    response = client.post("/models/wrong_id/simulate", json={'message': {}})
    assert response.status_code == 404


def test_simulate_no_operations(client):
    response = client.post("/models/iJO1366/simulate", json={})
    assert response.status_code == 200


def test_simulate_infeasible(client):
    measurements = [{'id': 'ATPM', 'measurements': [100, 100], 'type': 'reaction', 'db_name': 'bigg.reaction'}]
    response = client.post("/deltas", json={'model_id': 'iJO1366', 'conditions': {'measurements': measurements}})
    assert response.status_code == 200

    operations = response.json['operations']
    response = client.post("/models/iJO1366/simulate", json={'operations': operations})
    assert response.status_code == 200
    assert response.json['flux_distribution']['ATPM'] == pytest.approx(100)


def test_simulate_modify(monkeypatch, client):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, 'get_reaction_equations', lambda self, genotype: {})

    # Mock id-mapper api queries for efficiency
    def query_identifiers(object_ids, db_from, db_to):
        q = (object_ids, db_from, db_to)
        if q == (['MNXM2029', 'MNXM3447', 'MNXM368', 'MNXM7019'], 'mnx', 'chebi'):
            return {'MNXM368': ['16929', '10647', '12842', '26699', '52330', '57952']}
        elif q == (['MNXM2029', 'MNXM3447', 'MNXM368', 'MNXM7019'], 'mnx', 'bigg'):
            return {'MNXM3447': ['2agpe141'], 'MNXM368': ['g3pe'], 'MNXM7019': ['apg141'], 'MNXM2029': ['pg141']}
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
            return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
            return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
        elif q == (['glc__D', 'caro'], 'bigg', 'mnx'):
            return {'glc__D': ['MNXM41'], 'caro': ['MNXM614']}
        elif q == (['glc__D', 'caro'], 'bigg', 'chebi'):
            return {'glc__D': ['17634', '12965', '20999', '4167'], 'caro': ['17579', '10355', '12392', '22834', '40987']}  # noqa
        raise NotImplementedError(f"Unmocked query!")
    monkeypatch.setattr(adapter, 'query_identifiers', query_identifiers)

    for query, message in {'modify': MESSAGE_MODIFY, 'fluxes': MESSAGE_FLUXES}.items():
        response = client.post("/models/iJO1366/simulate", json={'message': message})
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


def test_simulate_different_objective(client):
    response = client.post("/models/iJO1366/simulate", json={'objective': 'EX_etoh_e'})
    assert response.status_code == 200
    result = response.json
    assert abs(result['flux_distribution']['EX_etoh_e']) - 20.0 < 0.001

    response = client.post("/models/iJO1366/simulate", json={})
    assert response.status_code == 200
    result = response.json
    assert abs(result['flux_distribution']['EX_etoh_e']) < 0.001


def test_deltas_post(monkeypatch, client):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, 'get_reaction_equations', lambda self, genotype: {})

    response = client.post('/deltas', json={
        'model_id': 'iJO1366',
        'conditions': {
            'medium': [
                {'id': 'chebi:44080'}, {'id': 'chebi:15075'}, {'id': 'chebi:15377'}, {'id': 'chebi:15378'}, {'id': 'chebi:15379'}, {'id': 'chebi:15982'}, {'id': 'chebi:16189'}, {'id': 'chebi:16526'}, {'id': 'chebi:16643'}, {'id': 'chebi:17883'}, {'id': 'chebi:18212'}, {'id': 'chebi:18367'}, {'id': 'chebi:18420'}, {'id': 'chebi:25371'}, {'id': 'chebi:27638'}, {'id': 'chebi:28938'}, {'id': 'chebi:29033'}, {'id': 'chebi:29034'}, {'id': 'chebi:29035'}, {'id': 'chebi:29036'}, {'id': 'chebi:29101'}, {'id': 'chebi:29103'}, {'id': 'chebi:29105'}, {'id': 'chebi:29108'}, {'id': 'chebi:36271'}, {'id': 'chebi:42758'}, {'id': 'chebi:49786'}  # noqa
            ],
            'genotype': ['+Aac', '-pta'],
            'measurements': [
                {'type': 'compound', 'id': 'chebi:42758', 'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'measurements': [-9.0]},
                {'type': 'compound', 'id': 'chebi:16236', 'unit': 'mmol', 'name': 'ethanol', 'measurements': [5.0, 4.8, 5.2, 4.9]},
                {'type': 'reaction', 'id': 'PFK', 'measurements': [5, 4.8, 7]},
                {'type': 'reaction', 'id': 'PGK', 'measurements': [5, 5]},
            ],
        },
        'operations': [{
            'operation': 'remove',
            'type': 'reaction',
            'id': 'KARA1',
        }],
    })
    assert response.status_code == 200
    assert len(response.json['operations']) == 72
