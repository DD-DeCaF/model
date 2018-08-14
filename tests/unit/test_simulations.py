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

from model.simulations import METHODS, simulate


@pytest.mark.skip(reason="TMY results is currently not implemented")
@pytest.mark.parametrize("objective", ['bigg:akg'])
def test_tmy_result(e_coli_core, objective):
    to_return = ['fluxes', 'tmy', 'model', 'growth-rate', 'removed-reactions']
    result = simulate(e_coli_core, 'fba', None, None, [objective], to_return)
    assert set(result) == set(to_return)


@pytest.mark.parametrize("method", METHODS)
def test_simulation_methods(e_coli_core, method):
    for method in METHODS:
        result = simulate(e_coli_core, method, None, None, [], None)
        if method not in {'fva', 'pfba-fva'}:
            reactions_ids = [i.id for i in e_coli_core.reactions]
            assert set(result['fluxes']) == set(reactions_ids)
