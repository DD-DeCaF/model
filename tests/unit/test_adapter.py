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
