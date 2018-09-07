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
from cobra.io.dict import model_to_dict

from model.ice_client import ICE


MEASUREMENTS = [
    {'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0], 'type': 'compound'},
    {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9], 'type': 'compound'},
    {'id': 'PFK', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'},
]


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
    response = client.post("/models/iJO1366/modify", json={'conditions': {'measurements': measurements}})
    assert response.status_code == 200

    operations = response.json['operations']
    response = client.post("/models/iJO1366/simulate", json={'operations': operations})
    assert response.status_code == 200
    assert response.json['flux_distribution']['ATPM'] == pytest.approx(100)


def test_simulate_fluxomics(monkeypatch, client):
    response = client.post("/models/iJO1366/modify", json={'conditions': {'measurements': MEASUREMENTS}})
    assert response.status_code == 200

    operations = response.json['operations']
    response = client.post("/models/iJO1366/simulate", json={'operations': operations})
    assert response.status_code == 200
    assert response.json['flux_distribution']['EX_glc__D_e'] == -9.0
    assert abs(response.json['flux_distribution']['EX_etoh_e'] - 4.64) < 0.001  # lower bound
    assert response.json['flux_distribution']['PFK'] == pytest.approx(5)


def test_simulate_modify(monkeypatch, client):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, 'get_reaction_equations', lambda self, genotype: {})

    conditions = {'measurements': MEASUREMENTS, 'genotype': ['+Aac', '-pta']}
    response = client.post("/models/iJO1366/modify", json={'conditions': conditions})
    assert response.status_code == 200

    operations = response.json['operations']
    assert any([op['operation'] == 'knockout' and op['type'] == 'gene' and op['id'] == 'b2297' for op in operations])
    assert any([op['operation'] == 'modify' and op['type'] == 'reaction' and op['id'] == 'EX_etoh_e' for op in operations])  # noqa
    assert any([op['operation'] == 'modify' and op['type'] == 'reaction' and op['id'] == 'PFK' for op in operations])

    response = client.post("/models/iJO1366/simulate", json={'method': 'pfba', 'operations': operations})
    assert response.status_code == 200
    fluxes = response.json['flux_distribution']

    assert fluxes['EX_glc__D_e'] == -9.0
    assert abs(fluxes['EX_etoh_e'] - 4.64) < 0.001  # lower bound
    assert fluxes['PFK'] == pytest.approx(5)


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

    response = client.post("/models/iJO1366/modify", json={
        'conditions': {
            'medium': [
                {'id': 'chebi:44080'}, {'id': 'chebi:15075'}, {'id': 'chebi:15377'}, {'id': 'chebi:15378'}, {'id': 'chebi:15379'}, {'id': 'chebi:15982'}, {'id': 'chebi:16189'}, {'id': 'chebi:16526'}, {'id': 'chebi:16643'}, {'id': 'chebi:17883'}, {'id': 'chebi:18212'}, {'id': 'chebi:18367'}, {'id': 'chebi:18420'}, {'id': 'chebi:25371'}, {'id': 'chebi:27638'}, {'id': 'chebi:28938'}, {'id': 'chebi:29033'}, {'id': 'chebi:29034'}, {'id': 'chebi:29035'}, {'id': 'chebi:29036'}, {'id': 'chebi:29101'}, {'id': 'chebi:29103'}, {'id': 'chebi:29105'}, {'id': 'chebi:29108'}, {'id': 'chebi:36271'}, {'id': 'chebi:42758'}, {'id': 'chebi:49786'}  # noqa
            ],
            'genotype': ['+Aac', '-pta'],
            'measurements': [
                {'type': 'compound', 'id': 'chebi:42758', 'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'measurements': [-9.0]},  # noqa
                {'type': 'compound', 'id': 'chebi:16236', 'unit': 'mmol', 'name': 'ethanol', 'measurements': [5.0, 4.8, 5.2, 4.9]},  # noqa
                {'type': 'reaction', 'id': 'PFK', 'measurements': [5, 4.8, 7]},
                {'type': 'reaction', 'id': 'PGK', 'measurements': [5, 5]},
            ],
        },
    })
    assert response.status_code == 200
    assert len(response.json['operations']) == 335


def test_simulate_custom_model(client, e_coli_core):
    model_serialized = model_to_dict(e_coli_core)
    response = client.post("/simulate", json={'model': model_serialized,
                                              'biomass_reaction': "BIOMASS_Ecoli_core_w_GAM"})
    assert response.status_code == 200
    assert response.json['growth_rate'] == pytest.approx(0.8739215069684307)
