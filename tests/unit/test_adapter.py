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

import pytest

from model.adapter import adapt_from_medium, adapt_from_genotype, adapt_from_measurements
from model.exceptions import NoIDMapping


def test_medium_adapter(iJO1366):
    medium = [
        {'id': 'chebi:63041'},
        {'id': 'chebi:91249'},
        {'id': 'chebi:86244'},
        {'id': 'chebi:131387'},
    ]
    operations = adapt_from_medium(iJO1366, medium)
    assert len(operations) == 38
    assert set(iJO1366.medium) == {'EX_fe3_e', 'EX_h2o_e', 'EX_mobd_e', 'EX_nh4_e', 'EX_so4_e', 'EX_ni2_e', 'EX_mn2_e', 'EX_cl_e'}  # noqa
    assert all(iJO1366.reactions.get_by_id(r).lower_bound == -1000 for r in iJO1366.medium)


def test_transport_reaction(iJO1366):
    changes = MeasurementChangeModel(iJO1366, [])
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
