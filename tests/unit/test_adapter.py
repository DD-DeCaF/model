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

import pytest

from model.adapter import (
    MeasurementChangeModel, MediumChangeModel, MediumSalts, NoIDMapping, get_unique_metabolite, next_measured_reaction)
from model.storage import Models


logging.disable(logging.CRITICAL)


@pytest.mark.asyncio
async def test_existing_metabolite():
    ecoli = Models.get('iJO1366')
    assert get_unique_metabolite(ecoli, 'chebi:17790') == get_unique_metabolite(
        ecoli, 'meoh', db_name='bigg.metabolite')
    assert get_unique_metabolite(ecoli, 'succ', db_name='bigg.metabolite').formula == 'C4H4O4'
    with pytest.raises(NoIDMapping):
        await get_unique_metabolite(ecoli, 'wrong_id')


def test_medium_salts():
    salts = MediumSalts.get()
    assert len(salts) > 2000
    assert len(salts['75832']) == 2
    assert len(salts['30808']) == 2
    assert len(salts['86254']) == 4


def test_medium_change_model():
    ecoli = Models.get('iJO1366')
    medium = [
        {'id': 'chebi:63041'},
        {'id': 'chebi:91249'},
        {'id': 'chebi:86244'},
        {'id': 'chebi:131387'},
    ]
    changes = MediumChangeModel(ecoli.copy(), medium)
    changes.apply_medium()
    model = changes.model
    assert 5 <= len(model.medium) <= 10
    assert {'EX_fe3_e', 'EX_h2o_e', 'EX_mobd_e', 'EX_nh4_e', 'EX_so4_e'} <= set(list(model.medium.keys()))


def test_transport_reaction():
    ecoli = Models.get('iJO1366')
    changes = MeasurementChangeModel(ecoli.copy(), [])
    assert changes.has_transport('o2', 1)
    assert changes.has_transport('fe2', -1)
    assert not changes.has_transport('btn', 1)
    changes.model.reactions.EX_btn_e.bounds = (0.1, 0.1)
    with pytest.warns(UserWarning):
        solution = changes.model.optimize()
    assert solution.status == 'infeasible'
    changes.allow_transport(changes.model.metabolites.btn_e, 1)
    assert changes.has_transport('btn', 1)
    solution = changes.model.optimize()
    assert solution.status == 'optimal'


def test_next_measured_reaction():
    ecoli = Models.get('iJO1366')
    assert next_measured_reaction(ecoli.reactions.EX_co2_e) == ecoli.reactions.CO2tex
    assert next_measured_reaction(ecoli.reactions.EX_glc__D_e) is None


def test_minimize_distance():
    model = Models.get('iJO1366')
    measurements = [{'type': 'reaction', 'id': 'GND', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.36, 2.45, 1.92]}, {'type': 'reaction', 'id': 'CS', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.5, 2.13, 1.54, 7.3]}, {'type': 'reaction', 'id': 'TPI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.31, 8.34, 8.4]}, {'type': 'reaction', 'id': 'FBA', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.31, 8.34, 8.4, 7.9]}, {'type': 'reaction', 'id': 'PFK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.31, 8.34, 8.4, 7.9]}, {'type': 'reaction', 'id': 'FUM', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.99, 1.59, 0.98, 6.7]}, {'type': 'reaction', 'id': 'GAPD', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [17.15, 17.13, 17.12, 16.8]}, {'type': 'reaction', 'id': 'G6PDH2r', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.55, 2.55, 2.08, 4.4]}, {'type': 'reaction', 'id': 'PGI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [7.3, 7.3, 7.75, 5.5]}, {'type': 'reaction', 'id': 'GLCptspp', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [10.0, 10.0, 10.0, 10.0]}, {'type': 'reaction', 'id': 'ICDHyr', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.5, 2.02, 0.86]}, {'type': 'reaction', 'id': 'ME1', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.26, 0.03, 0.03, 0.4]}, {'type': 'reaction', 'id': 'ME2', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.26, 0.03, 0.03, 0.4]}, {'type': 'reaction', 'id': 'MDH', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.42, 1.64, 1.04]}, {'type': 'reaction', 'id': 'PYK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.71, 3.02, 2.89, 12.3]}, {'type': 'reaction', 'id': 'PPC', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.63, 2.36, 2.56, 2.8]}, {'type': 'reaction', 'id': 'PPCK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.63, 2.36, 2.56, 2.8]}, {'type': 'reaction', 'id': 'PDHm', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [11.57, 11.2, 11.57, 9.4]}, {'type': 'reaction', 'id': 'TKT1', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.05, 1.09, 0.71, 1.4]}, {'type': 'reaction', 'id': 'RPI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.3, 1.36, 1.21]}, {'type': 'reaction', 'id': 'RBP4E', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.05, 1.09, 0.71]}, {'type': 'reaction', 'id': 'SUCDi', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.72, 1.32, 0.68]}, {'type': 'reaction', 'id': 'MALS', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.11, 0.68]}, {'type': 'reaction', 'id': 'ICL', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.11, 0.68]}, {'type': 'reaction', 'id': 'TALA', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.4]}, {'type': 'reaction', 'id': 'TKT2', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.1]}]  # noqa
    minimized_measurements = [{'type': 'reaction', 'id': 'GND', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.2433333333174232]}, {'type': 'reaction', 'id': 'CS', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [3.367499999993015]}, {'type': 'reaction', 'id': 'TPI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.34999999999418]}, {'type': 'reaction', 'id': 'FBA', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.237500000004367]}, {'type': 'reaction', 'id': 'PFK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.23750000000291]}, {'type': 'reaction', 'id': 'FUM', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.8150000000034923]}, {'type': 'reaction', 'id': 'GAPD', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [17.04999999998981]}, {'type': 'reaction', 'id': 'G6PDH2r', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.6033333333523387]}, {'type': 'reaction', 'id': 'PGI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [6.962500000020441]}, {'type': 'reaction', 'id': 'GLCptspp', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [10.0]}, {'type': 'reaction', 'id': 'ICDHyr', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.7933333333348855]}, {'type': 'reaction', 'id': 'ME1', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.17999999999301508]}, {'type': 'reaction', 'id': 'ME2', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.17999999999301508]}, {'type': 'reaction', 'id': 'MDH', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [-1.3666666666666742]}, {'type': 'reaction', 'id': 'PYK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [5.2299999999813735]}, {'type': 'reaction', 'id': 'PPC', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.5874999999941792]}, {'type': 'reaction', 'id': 'PPCK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.5874999999941792]}, {'type': 'reaction', 'id': 'PDHm', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [11.57, 11.2, 11.57, 9.4]}, {'type': 'reaction', 'id': 'TKT1', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.0625]}, {'type': 'reaction', 'id': 'RPI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [-1.2899999999708962]}, {'type': 'reaction', 'id': 'RBP4E', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.0]}, {'type': 'reaction', 'id': 'SUCDi', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.2399999999906868]}, {'type': 'reaction', 'id': 'MALS', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.9133333333303697]}, {'type': 'reaction', 'id': 'ICL', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.0]}, {'type': 'reaction', 'id': 'TALA', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.0625]}, {'type': 'reaction', 'id': 'TKT2', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [-0.10916666665347297]}]
    changes = MeasurementChangeModel(model, measurements)
    changes.minimize_distance()
    for i, change in enumerate(changes.measurements):
        assert change['measurements'] == pytest.approx(minimized_measurements[i]['measurements'])
