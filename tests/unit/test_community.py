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

from simulations.modeling.community import METHODS, simulate


@pytest.mark.slow
@pytest.mark.parametrize("method", METHODS)
def test_community_simulation(e_coli_core, iJO1366, method):
    e_coli_core, biomass_reaction, is_ec_model = e_coli_core
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    result = simulate([e_coli_core, iJO1366], [], method)
    assert result["growth_rate"] == pytest.approx(0)
    assert result["abundance"]["e_coli_core"] == pytest.approx(0.967027244770635)
    assert result["abundance"]["iJO1366"] == pytest.approx(0.03297275522936505)
