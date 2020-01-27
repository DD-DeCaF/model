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

from simulations.modeling.pytfa_helpers import HandlerThermo


@pytest.fixture(scope="function")
def thermo_model(iJO1366):
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    tmodel = HandlerThermo(iJO1366)
    tmodel._convert()
    return tmodel


def test_handlerThermo_consistency(iJO1366):
    """Check if FBA isn't altered after the conversion."""
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    fba_solution = iJO1366.optimize().fluxes
    tmodel = HandlerThermo(iJO1366)
    tmodel._convert()
    assert (tmodel.optimize().fluxes == fba_solution).all()


def test_tmfa(iJO1366, tmodel):
    """Check if it affects the solution."""
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    fba_solution = iJO1366.optimize().fluxes
    thermo_solution = tmodel.tmfa().fluxes
    assert (thermo_solution != fba_solution).any()


def test_metabolomics(tmodel):
    """Check if it affects the TMFA solution."""
    solution_before = tmodel.tmfa().fluxes
    metabolomics = [
        {
            "name": "D-Glucose",
            "identifier": "glc__D",
            "namespace": "bigg.metabolite",
            "measurement": 0.2,
            "uncertainty": 0.01,
        }
    ]
    tmodel._apply_metabolomics(metabolomics)
    solution_metabolomics = tmodel.tmfa().fluxes
    assert (solution_before != solution_metabolomics).any()
