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

"""Test local HTTP endpoints."""

from collections import namedtuple

import pytest
import requests

from simulations.ice_client import ICE


FLUXOMICS = [
    {
        "name": "Phosphofructokinase",
        "identifier": "PFK",
        "namespace": "bigg.reaction",
        "measurement": 5,
        "uncertainty": 0,
    }
]

UPTAKE_SECRETION_RATES = [
    {
        "name": "aldehydo-D-glucose",
        "identifier": "CHEBI:42758",
        "namespace": "chebi",
        "measurement": -9.0,
        "uncertainty": 0,
    },
    {
        "name": "ethanol",
        "identifier": "CHEBI:16236",
        "namespace": "chebi",
        "measurement": 4.9,
        "uncertainty": 0,
    },
]


def test_simulate_wrong_id(monkeypatch, client):
    # Mock `requests` to skip the external API request
    Response = namedtuple("Response", ["status_code"])
    monkeypatch.setattr(
        requests, "get", lambda *args, **kwargs: Response(status_code=404)
    )
    response = client.post("/simulate", json={"model_id": 404, "message": {}})
    assert response.status_code == 404


def test_simulate_unauthorized(client, models):
    response = client.post(
        "/simulate", json={"model_id": models["e_coli_core_proprietary"]}
    )
    assert response.status_code == 403


def test_simulate_no_operations(client, models):
    response = client.post("/simulate", json={"model_id": models["iJO1366"]})
    assert response.status_code == 200
    assert response.json["status"] == "optimal"


def test_simulate_infeasible(client, models):
    fluxomics = [
        {
            "name": "E. coli biomass objective function",
            "identifier": "BIOMASS_Ec_iJO1366_core_53p95M",
            "namespace": "bigg.reaction",
            # Force an impossible growth to ensure infeasability
            "measurement": 1000,
            "uncertainty": 0,
        }
    ]
    response = client.post(
        f"/models/{models['iJO1366']}/modify", json={"fluxomics": fluxomics}
    )
    assert response.status_code == 200

    operations = response.json["operations"]
    response = client.post(
        "/simulate", json={"model_id": models["iJO1366"], "operations": operations}
    )
    assert response.status_code == 200
    assert response.json["status"] == "infeasible"


def test_simulate_fluxomics(monkeypatch, client, models):
    response = client.post(
        f"/models/{models['iJO1366']}/modify", json={"fluxomics": FLUXOMICS}
    )
    assert response.status_code == 200

    operations = response.json["operations"]
    response = client.post(
        "/simulate", json={"model_id": models["iJO1366"], "operations": operations}
    )
    assert response.status_code == 200
    assert response.json["status"] == "optimal"


def test_simulate_modify(monkeypatch, client, models):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, "get_reaction_equations", lambda self, genotype: {})

    experiment = {
        "fluxomics": FLUXOMICS,
        "uptake_secretion_rates": UPTAKE_SECRETION_RATES,
        "genotype": "+Aac,-pta",
    }
    response = client.post(f"/models/{models['iJO1366']}/modify", json=experiment)
    assert response.status_code == 200

    operations = response.json["operations"]
    assert any(
        [
            op["operation"] == "knockout"
            and op["type"] == "gene"
            and op["id"] == "b2297"
            for op in operations
        ]
    )
    assert any(
        [
            op["operation"] == "modify"
            and op["type"] == "reaction"
            and op["id"] == "EX_etoh_e"
            for op in operations
        ]
    )
    assert any(
        [
            op["operation"] == "modify"
            and op["type"] == "reaction"
            and op["id"] == "PFK"
            for op in operations
        ]
    )

    response = client.post(
        "/simulate",
        json={
            "model_id": models["iJO1366"],
            "method": "pfba",
            "operations": operations,
        },
    )
    assert response.status_code == 200
    assert response.json["status"] == "optimal"
    fluxes = response.json["flux_distribution"]

    assert fluxes["EX_glc__D_e"] == -9.0
    assert fluxes["PFK"] == pytest.approx(5)


def test_simulate_different_objective(client, models):
    response = client.post(
        "/simulate", json={"model_id": models["iJO1366"], "objective_id": "EX_etoh_e"}
    )
    assert response.status_code == 200
    result = response.json
    assert result["status"] == "optimal"
    assert abs(result["flux_distribution"]["EX_etoh_e"]) == pytest.approx(20)

    response = client.post("/simulate", json={"model_id": models["iJO1366"]})
    assert response.status_code == 200
    result = response.json
    assert result["status"] == "optimal"
    assert abs(result["flux_distribution"]["EX_etoh_e"]) == pytest.approx(0)


def test_modify(monkeypatch, client, models):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, "get_reaction_equations", lambda self, genotype: {})

    response = client.post(
        f"/models/{models['iJO1366']}/modify",
        json={
            "medium": [
                {
                    "name": "methanol",
                    "identifier": "CHEBI:44080",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "selenate",
                    "identifier": "CHEBI:15075",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "water",
                    "identifier": "CHEBI:15377",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "hydron",
                    "identifier": "CHEBI:15378",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "dioxygen",
                    "identifier": "CHEBI:15379",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "cob(I)alamin",
                    "identifier": "CHEBI:15982",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "sulfate",
                    "identifier": "CHEBI:16189",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "carbon dioxide",
                    "identifier": "CHEBI:16526",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "L-methionine",
                    "identifier": "CHEBI:16643",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "hydrogen chloride",
                    "identifier": "CHEBI:17883",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "selenite(2-)",
                    "identifier": "CHEBI:18212",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "phosphate(3-)",
                    "identifier": "CHEBI:18367",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "magnesium(2+)",
                    "identifier": "CHEBI:18420",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "molybdic acid",
                    "identifier": "CHEBI:25371",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "cobalt atom",
                    "identifier": "CHEBI:27638",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "ammonium",
                    "identifier": "CHEBI:28938",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "iron(2+)",
                    "identifier": "CHEBI:29033",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "iron(3+)",
                    "identifier": "CHEBI:29034",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "manganese(2+)",
                    "identifier": "CHEBI:29035",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "copper(2+)",
                    "identifier": "CHEBI:29036",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "sodium(1+)",
                    "identifier": "CHEBI:29101",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "potassium(1+)",
                    "identifier": "CHEBI:29103",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "zinc(2+)",
                    "identifier": "CHEBI:29105",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "calcium(2+)",
                    "identifier": "CHEBI:29108",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "hydrogentungstate",
                    "identifier": "CHEBI:36271",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "aldehydo-D-glucose",
                    "identifier": "CHEBI:42758",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
                {
                    "name": "nickel(2+)",
                    "identifier": "CHEBI:49786",
                    "namespace": "chebi",
                    "mass_concentration": None,
                },
            ],
            "genotype": "+Aac,-pta",
            "uptake_secretion_rates": [
                {
                    "name": "aldehydo-D-glucose",
                    "identifier": "CHEBI:42758",
                    "namespace": "chebi",
                    "measurement": -9.0,
                    "uncertainty": 0,
                },
                {
                    "name": "ethanol",
                    "identifier": "CHEBI:16236",
                    "namespace": "chebi",
                    "measurement": 4.9,
                    "uncertainty": 0,
                },
            ],
            "fluxomics": [
                {
                    "name": "Phosphofructokinase",
                    "identifier": "PFK",
                    "namespace": "bigg.reaction",
                    "measurement": 4.8,
                    "uncertainty": 0,
                },
                {
                    "name": "Phosphoglycerate kinase",
                    "identifier": "PGK",
                    "namespace": "bigg.reaction",
                    "measurement": 5,
                    "uncertainty": 0,
                },
            ],
        },
    )
    assert response.status_code == 200
    assert len(response.json["operations"]) == 329


def test_prokaryomics_md120_bw25113(client, models):
    """Test constraining and simulating a model with a real data set."""
    data = {
        "medium": [
            {
                "name": "dipotassium hydrogen phosphate",
                "identifier": "CHEBI:131527",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "calcium dichloride",
                "identifier": "CHEBI:3312",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "L-glutamic acid",
                "identifier": "CHEBI:16015",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "sodium chloride",
                "identifier": "CHEBI:26710",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "iron trichloride",
                "identifier": "CHEBI:30808",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "magnesium sulfate",
                "identifier": "CHEBI:32599",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "disodium hydrogenphosphate",
                "identifier": "CHEBI:34683",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "cobalt dichloride",
                "identifier": "CHEBI:35696",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "aldehydo-D-glucose",
                "identifier": "CHEBI:42758",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "copper(II) chloride",
                "identifier": "CHEBI:49553",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "zinc dichloride",
                "identifier": "CHEBI:49976",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "ammonium sulfate",
                "identifier": "CHEBI:62946",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "manganese(II) chloride",
                "identifier": "CHEBI:63041",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "potassium nitrate",
                "identifier": "CHEBI:63043",
                "namespace": "chebi",
                "mass_concentration": None,
            },
            {
                "name": "sodium molybdate (anhydrous)",
                "identifier": "CHEBI:75215",
                "namespace": "chebi",
                "mass_concentration": None,
            },
        ],
        "fluxomics": [
            {
                "name": "Citrate synthase",
                "identifier": "CS",
                "namespace": "bigg.reaction",
                "measurement": 7.2,
                "uncertainty": 0,
            },
            {
                "name": "Fructose-bisphosphate aldolase",
                "identifier": "FBA",
                "namespace": "bigg.reaction",
                "measurement": 7.9,
                "uncertainty": 0,
            },
            {
                "name": "Phosphofructokinase",
                "identifier": "PFK",
                "namespace": "bigg.reaction",
                "measurement": 7.9,
                "uncertainty": 0,
            },
            {
                "name": "Fumarase",
                "identifier": "FUM",
                "namespace": "bigg.reaction",
                "measurement": 6.7,
                "uncertainty": 0,
            },
            {
                "name": "Glyceraldehyde-3-phosphate dehydrogenase",
                "identifier": "GAPD",
                "namespace": "bigg.reaction",
                "measurement": 16.8,
                "uncertainty": 0,
            },
            {
                "name": "Glucose 6-phosphate dehydrogenase",
                "identifier": "G6PDH2r",
                "namespace": "bigg.reaction",
                "measurement": 4.6,
                "uncertainty": 0,
            },
            {
                "name": "Glucose-6-phosphate isomerase",
                "identifier": "PGI",
                "namespace": "bigg.reaction",
                "measurement": 5.3,
                "uncertainty": 0,
            },
            {
                "name": "D-glucose transport via PEP:Pyr PTS (periplasm)",
                "identifier": "GLCptspp",
                "namespace": "bigg.reaction",
                "measurement": 10.0,
                "uncertainty": 0,
            },
            {
                "name": "Malic enzyme (NAD)",
                "identifier": "ME1",
                "namespace": "bigg.reaction",
                "measurement": 0.3,
                "uncertainty": 0,
            },
            {
                "name": "Malic enzyme (NADP)",
                "identifier": "ME2",
                "namespace": "bigg.reaction",
                "measurement": 0.3,
                "uncertainty": 0,
            },
            {
                "name": "Malate dehydrogenase",
                "identifier": "MDH",
                "namespace": "bigg.reaction",
                "measurement": 6.5,
                "uncertainty": 0,
            },
            {
                "name": "Pyruvate kinase",
                "identifier": "PYK",
                "namespace": "bigg.reaction",
                "measurement": 12.0,
                "uncertainty": 0,
            },
            {
                "name": "Phosphoenolpyruvate carboxylase",
                "identifier": "PPC",
                "namespace": "bigg.reaction",
                "measurement": 3.0,
                "uncertainty": 0,
            },
            {
                "name": "Phosphoenolpyruvate carboxykinase",
                "identifier": "PPCK",
                "namespace": "bigg.reaction",
                "measurement": 3.0,
                "uncertainty": 0,
            },
            {
                "name": "Pyruvate dehydrogenase",
                "identifier": "PDH",
                "namespace": "bigg.reaction",
                "measurement": 9.4,
                "uncertainty": 0,
            },
            {
                "name": "Transketolase",
                "identifier": "TKT1",
                "namespace": "bigg.reaction",
                "measurement": 1.5,
                "uncertainty": 0,
            },
            {
                "name": "Transaldolase",
                "identifier": "TALA",
                "namespace": "bigg.reaction",
                "measurement": 1.5,
                "uncertainty": 0,
            },
            {
                "name": "Transketolase",
                "identifier": "TKT2",
                "namespace": "bigg.reaction",
                "measurement": 1.1,
                "uncertainty": 0,
            },
        ],
        "genotype": "-b3643,-b0062,-b0063,-b0061,-b4350,-b3902,-b3903",
    }

    response = client.post(f"/models/{models['iJO1366']}/modify", json=data)
    assert response.status_code == 200

    response = client.post(
        "/simulate",
        json={"model_id": models["iJO1366"], "operations": response.json["operations"]},
    )
    assert response.status_code == 200
    assert response.json["status"] == "optimal"
    assert response.json["growth_rate"] == pytest.approx(0.5134445454218568)


def test_growth_rate_measurement(client, models):
    """Constrain the model with a single growth rate measurement."""
    data = {"growth_rate": {"measurement": 0.3, "uncertainty": 0}}

    response = client.post(f"/models/{models['iJO1366']}/modify", json=data)
    assert response.status_code == 200

    response = client.post(
        "/simulate",
        json={"model_id": models["iJO1366"], "operations": response.json["operations"]},
    )
    assert response.status_code == 200
    assert response.json["status"] == "optimal"
    assert response.json["growth_rate"] == pytest.approx(0.3)
