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

import gnomic

from simulations.ice_client import ICE
from simulations.modeling.adapter import (
    SALTS,
    apply_genotype,
    apply_measurements,
    apply_medium,
)


def test_medium_salts():
    assert len(SALTS) > 2000
    assert len(SALTS["CHEBI:75832"]["ions"]) == 2
    assert len(SALTS["CHEBI:30808"]["metals"]) == 7
    assert len(SALTS["CHEBI:86254"]["ions"]) == 2


def test_medium_adapter(iJO1366):
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    medium = [
        {"name": "Foo", "identifier": "CHEBI:63041", "namespace": "chebi"},
        {"name": "Bar", "identifier": "CHEBI:91249", "namespace": "chebi"},
        {"name": "Baz", "identifier": "CHEBI:86244", "namespace": "chebi"},
        {"name": "Goo", "identifier": "CHEBI:131387", "namespace": "chebi"},
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
    }
    assert all(
        iJO1366.reactions.get_by_id(r).lower_bound == -1000 for r in iJO1366.medium
    )


def test_genotype_adapter(monkeypatch, iJO1366):
    iJO1366, biomass_reaction, is_ec_model = iJO1366

    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, "get_reaction_equations", lambda self, genotype: {})

    genotype_changes = gnomic.Genotype.parse("+Aac,-pta")
    operations, warnings, errors = apply_genotype(iJO1366, genotype_changes)
    assert len(operations) == 1
    assert len(errors) == 0


def test_measurements_adapter(iJO1366):
    iJO1366, biomass_reaction, is_ec_model = iJO1366
    uptake_secretion_rates = [
        {
            "name": "Foo",
            "identifier": "CHEBI:42758",
            "namespace": "chebi",
            "measurement": -9.0,
            "uncertainty": 0,
        },
        {
            "name": "Foo",
            "identifier": "CHEBI:16236",
            "namespace": "chebi",
            "measurement": 4.9,
            "uncertainty": 0,
        },
    ]
    fluxomics = [
        {
            "name": "Foo",
            "identifier": "PFK",
            "namespace": "bigg.reaction",
            "measurement": 4.8,
            "uncertainty": 0,
        },
        {
            "name": "Foo",
            "identifier": "PGK",
            "namespace": "bigg.reaction",
            "measurement": 5,
            "uncertainty": 0,
        },
    ]
    proteomics = {"identifier": "P0A8V2", "measurement": 5.03e-6, "uncertainty": 0}
    operations, warnings, errors = apply_measurements(
        iJO1366,
        biomass_reaction,
        is_ec_model,
        fluxomics,
        [],
        proteomics,
        uptake_secretion_rates,
        [],
        None,
    )
    # 4 operations (2 rates + 2 fluxomics) + 1 warning (not an ecModel) are expected:
    assert len(operations) == 4
    assert len(warnings) == 1
    assert len(errors) == 0


def test_proteomics_adapter(eciML1515):
    eciML1515, biomass_reaction, is_ec_model = eciML1515
    proteomics = [
        {
            "identifier": "P0A8V2",  # protein not in model (should be skipped)
            "measurement": 5.03e-6,
            "uncertainty": 0,
        },
        {
            "identifier": "P0AFG8",
            "measurement": 8.2e-3,  # very high value (should be kept)
            "uncertainty": 8.2e-6,
        },
        {
            "identifier": "P15254",
            "measurement": 6.54e-8,  # very low value (should be removed)
            "uncertainty": 0,
        },
        {
            "identifier": "P0A6C5",
            "measurement": 5.93e-8,  # very low value (should be removed)
            "uncertainty": 0,
        },
    ]
    growth_rate = {"measurement": 0.1, "uncertainty": 0.01}
    operations, warnings, errors = apply_measurements(
        eciML1515,
        biomass_reaction,
        is_ec_model,
        [],
        [],
        proteomics,
        [],
        [],
        growth_rate,
    )
    # 2 operations (1 proteomics + growth rate) + 1 warning (1 skipped protein) are
    # expected:
    assert len(operations) == 2
    assert len(warnings) == 1
    assert len(errors) == 0
