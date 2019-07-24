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

import numpy as np
import pandas as pd
import pytest

from simulations.modeling.driven import adjust_fluxes2model, minimize_distance


def test_minimize_distance_no_growth_rate(iJO1366):
    iJO1366, biomass_reaction = iJO1366
    with pytest.raises(ValueError):
        measurements = [
            {
                "type": "reaction",
                "id": "GND",
                "namespace": "bigg.reaction",
                "measurements": [2.36, 2.45, 1.92],
            }
        ]
        minimize_distance(iJO1366, biomass_reaction, None, measurements)


def test_adjust_fluxes2model(iJO1366):
    iJO1366, biomass_reaction = iJO1366

    # Since there is not a general flux direction of the system imposed (via
    # biomass or other constraint), we artificially increase the lb of the
    # biomass reaction in order to force the flux direction.
    # Consider replacing the measurements with an experiment that does include a
    # growth rate constraint.
    iJO1366.reactions.BIOMASS_Ec_iJO1366_WT_53p95M.lower_bound = 0.5

    measurements = [
        {
            "type": "reaction",
            "id": "GND",
            "namespace": "bigg.reaction",
            "measurements": [2.36, 2.45, 1.92],
        },
        {
            "type": "reaction",
            "id": "CS",
            "namespace": "bigg.reaction",
            "measurements": [2.5, 2.13, 1.54, 7.3],
        },
        {
            "type": "reaction",
            "id": "TPI",
            "namespace": "bigg.reaction",
            "measurements": [8.31, 8.34, 8.4],
        },
        {
            "type": "reaction",
            "id": "FBA",
            "namespace": "bigg.reaction",
            "measurements": [8.31, 8.34, 8.4, 7.9],
        },
        {
            "type": "reaction",
            "id": "PFK",
            "namespace": "bigg.reaction",
            "measurements": [8.31, 8.34, 8.4, 7.9],
        },
        {
            "type": "reaction",
            "id": "FUM",
            "namespace": "bigg.reaction",
            "measurements": [1.99, 1.59, 0.98, 6.7],
        },
        {
            "type": "reaction",
            "id": "GAPD",
            "namespace": "bigg.reaction",
            "measurements": [17.15, 17.13, 17.12, 16.8],
        },
        {
            "type": "reaction",
            "id": "G6PDH2r",
            "namespace": "bigg.reaction",
            "measurements": [2.55, 2.55, 2.08, 4.4],
        },
        {
            "type": "reaction",
            "id": "PGI",
            "namespace": "bigg.reaction",
            "measurements": [7.3, 7.3, 7.75, 5.5],
        },
        {
            "type": "reaction",
            "id": "GLCptspp",
            "namespace": "bigg.reaction",
            "measurements": [10.0, 10.0, 10.0, 10.0],
        },
        {
            "type": "reaction",
            "id": "ICDHyr",
            "namespace": "bigg.reaction",
            "measurements": [2.5, 2.02, 0.86],
        },
        {
            "type": "reaction",
            "id": "ME1",
            "namespace": "bigg.reaction",
            "measurements": [0.26, 0.03, 0.03, 0.4],
        },
        {
            "type": "reaction",
            "id": "ME2",
            "namespace": "bigg.reaction",
            "measurements": [0.26, 0.03, 0.03, 0.4],
        },
        {
            "type": "reaction",
            "id": "MDH",
            "namespace": "bigg.reaction",
            "measurements": [1.42, 1.64, 1.04],
        },
        {
            "type": "reaction",
            "id": "PYK",
            "namespace": "bigg.reaction",
            "measurements": [2.71, 3.02, 2.89, 12.3],
        },
        {
            "type": "reaction",
            "id": "PPC",
            "namespace": "bigg.reaction",
            "measurements": [2.63, 2.36, 2.56, 2.8],
        },
        {
            "type": "reaction",
            "id": "PPCK",
            "namespace": "bigg.reaction",
            "measurements": [2.63, 2.36, 2.56, 2.8],
        },
        {
            "type": "reaction",
            "id": "PDH",
            "namespace": "bigg.reaction",
            "measurements": [11.57, 11.2, 11.57, 9.4],
        },
        {
            "type": "reaction",
            "id": "TKT1",
            "namespace": "bigg.reaction",
            "measurements": [1.05, 1.09, 0.71, 1.4],
        },
        {
            "type": "reaction",
            "id": "RPI",
            "namespace": "bigg.reaction",
            "measurements": [1.3, 1.36, 1.21],
        },
        {
            "type": "reaction",
            "id": "RBP4E",
            "namespace": "bigg.reaction",
            "measurements": [1.05, 1.09, 0.71],
        },
        {
            "type": "reaction",
            "id": "SUCDi",
            "namespace": "bigg.reaction",
            "measurements": [1.72, 1.32, 0.68],
        },
        {
            "type": "reaction",
            "id": "MALS",
            "namespace": "bigg.reaction",
            "measurements": [0.11, 0.68],
        },
        {
            "type": "reaction",
            "id": "ICL",
            "namespace": "bigg.reaction",
            "measurements": [0.11, 0.68],
        },
        {
            "type": "reaction",
            "id": "TALA",
            "namespace": "bigg.reaction",
            "measurements": [1.4],
        },
        {
            "type": "reaction",
            "id": "TKT2",
            "namespace": "bigg.reaction",
            "measurements": [1.1],
        },
    ]
    expected_distance = 4.269723695500828

    index = []
    observations = []
    for measure in measurements:
        index.append(measure["id"])
        observations.append(np.mean(measure["measurements"]))
    observations = pd.Series(index=index, data=observations)
    solution = adjust_fluxes2model(iJO1366, observations)
    minimized_fluxes = solution.fluxes.to_dict()

    # Calculate the total distance from the observed fluxes
    measurements = {d["id"]: np.mean(d["measurements"]) for d in measurements}
    calculated_distance = sum(
        [
            abs(abs(minimized_fluxes[reaction_id]) - measured_flux)
            for reaction_id, measured_flux in measurements.items()
        ]
    )

    assert expected_distance == pytest.approx(solution.objective_value)
    assert expected_distance == pytest.approx(calculated_distance)
