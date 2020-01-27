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
from simulations.storage import ModelWrapper


@pytest.mark.slow
@pytest.mark.parametrize("method", METHODS)
def test_community_simulation(e_coli_core, iJO1366, method):
    e_coli_core, biomass_reaction, is_ec_model = e_coli_core
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    # Create fake wrappers
    wrappers = [
        ModelWrapper(
            id=1,
            model=e_coli_core,
            project_id=None,
            organism_id=None,
            biomass_reaction=None,
            is_ec_model=False,
        ),
        ModelWrapper(
            id=2,
            model=iJO1366,
            project_id=None,
            organism_id=None,
            biomass_reaction=None,
            is_ec_model=False,
        ),
    ]
    result = simulate(wrappers, [], method)
    assert result["growth_rate"] == pytest.approx(0)
    abundance = next(r for r in result["abundance"] if r["id"] == "e_coli_core")
    assert abundance["value"] == pytest.approx(0.967027244770635)
    abundance = next(r for r in result["abundance"] if r["id"] == "iJO1366")
    assert abundance["value"] == pytest.approx(0.03297275522936505)
