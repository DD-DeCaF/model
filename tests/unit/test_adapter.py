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

from simulations.ice_client import ICE
from simulations.modeling.adapter import (
    SALTS,
    apply_genotype,
    apply_measurements,
    apply_medium,
)
from simulations.modeling.gnomic_helpers import full_genotype


def test_medium_salts():
    assert len(SALTS) > 2000
    assert len(SALTS["CHEBI:75832"]["ions"]) == 2
    assert len(SALTS["CHEBI:30808"]["metals"]) == 7
    assert len(SALTS["CHEBI:86254"]["ions"]) == 2


def test_medium_adapter(iJO1366):
    iJO1366, biomass_reaction = iJO1366
    medium = [
        {"id": "CHEBI:63041", "namespace": "chebi"},
        {"id": "CHEBI:91249", "namespace": "chebi"},
        {"id": "CHEBI:86244", "namespace": "chebi"},
        {"id": "CHEBI:131387", "namespace": "chebi"},
    ]
    operations, warnings, errors = apply_medium(iJO1366, medium)
    # 30 warnings are expected; 29 unique compounds not found in the model and 1
    # unmapped ion from the salts mapping
    assert len(warnings) == 30
    assert len(errors) == 0
    assert set(iJO1366.medium) == {
        "EX_fe2_e",
        "EX_fe3_e",
        "EX_h2o_e",
        "EX_mobd_e",
        "EX_nh4_e",
        "EX_so4_e",
        "EX_ni2_e",
        "EX_mn2_e",
        "EX_cl_e",
    }  # noqa
    assert all(
        iJO1366.reactions.get_by_id(r).lower_bound == -1000 for r in iJO1366.medium
    )


def test_genotype_adapter(monkeypatch, iJO1366):
    iJO1366, biomass_reaction = iJO1366

    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, "get_reaction_equations", lambda self, genotype: {})

    genotype_changes = full_genotype(["+Aac", "-pta"])
    operations, warnings, errors = apply_genotype(iJO1366, genotype_changes)
    assert len(operations) == 1
    assert len(errors) == 0


def test_measurements_adapter(iJO1366):
    iJO1366, biomass_reaction = iJO1366
    measurements = [
        {
            "type": "compound",
            "id": "CHEBI:42758",
            "namespace": "chebi",
            "measurements": [-9.0],
        },
        {
            "type": "compound",
            "id": "CHEBI:16236",
            "namespace": "chebi",
            "measurements": [5.0, 4.8, 5.2, 4.9],
        },
        {
            "type": "reaction",
            "id": "PFK",
            "namespace": "bigg.reaction",
            "measurements": [5, 4.8, 7],
        },
        {
            "type": "reaction",
            "id": "PGK",
            "namespace": "bigg.reaction",
            "measurements": [5, 5],
        },
    ]
    operations, warnings, errors = apply_measurements(
        iJO1366, biomass_reaction, None, measurements
    )
    assert len(operations) == 4
    assert len(errors) == 0
