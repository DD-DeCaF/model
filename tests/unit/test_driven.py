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

from model.driven import adjust_fluxes2model


def test_driven(iJO1366):
    measurements = [{'type': 'reaction', 'id': 'GND', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.36, 2.45, 1.92]}, {'type': 'reaction', 'id': 'CS', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.5, 2.13, 1.54, 7.3]}, {'type': 'reaction', 'id': 'TPI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.31, 8.34, 8.4]}, {'type': 'reaction', 'id': 'FBA', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.31, 8.34, 8.4, 7.9]}, {'type': 'reaction', 'id': 'PFK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [8.31, 8.34, 8.4, 7.9]}, {'type': 'reaction', 'id': 'FUM', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.99, 1.59, 0.98, 6.7]}, {'type': 'reaction', 'id': 'GAPD', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [17.15, 17.13, 17.12, 16.8]}, {'type': 'reaction', 'id': 'G6PDH2r', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.55, 2.55, 2.08, 4.4]}, {'type': 'reaction', 'id': 'PGI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [7.3, 7.3, 7.75, 5.5]}, {'type': 'reaction', 'id': 'GLCptspp', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [10.0, 10.0, 10.0, 10.0]}, {'type': 'reaction', 'id': 'ICDHyr', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.5, 2.02, 0.86]}, {'type': 'reaction', 'id': 'ME1', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.26, 0.03, 0.03, 0.4]}, {'type': 'reaction', 'id': 'ME2', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.26, 0.03, 0.03, 0.4]}, {'type': 'reaction', 'id': 'MDH', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.42, 1.64, 1.04]}, {'type': 'reaction', 'id': 'PYK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.71, 3.02, 2.89, 12.3]}, {'type': 'reaction', 'id': 'PPC', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.63, 2.36, 2.56, 2.8]}, {'type': 'reaction', 'id': 'PPCK', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [2.63, 2.36, 2.56, 2.8]}, {'type': 'reaction', 'id': 'PDHm', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [11.57, 11.2, 11.57, 9.4]}, {'type': 'reaction', 'id': 'TKT1', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.05, 1.09, 0.71, 1.4]}, {'type': 'reaction', 'id': 'RPI', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.3, 1.36, 1.21]}, {'type': 'reaction', 'id': 'RBP4E', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.05, 1.09, 0.71]}, {'type': 'reaction', 'id': 'SUCDi', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.72, 1.32, 0.68]}, {'type': 'reaction', 'id': 'MALS', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.11, 0.68]}, {'type': 'reaction', 'id': 'ICL', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [0.11, 0.68]}, {'type': 'reaction', 'id': 'TALA', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.4]}, {'type': 'reaction', 'id': 'TKT2', 'db_name': 'bigg.reaction', 'mode': 'quantitative', 'measurements': [1.1]}]  # noqa
    expected_distance = 3.1458333332930755

    index = []
    observations = []
    for measure in measurements:
        index.append(measure['id'])
        observations.append(np.mean(measure['measurements']))
    observations = pd.Series(index=index, data=observations)
    solution = adjust_fluxes2model(iJO1366, observations)
    minimized_fluxes = solution.fluxes.to_dict()

    # Calculate the total distance from the observed fluxes
    measurements = {d['id']: np.mean(d['measurements']) for d in measurements}
    calculated_distance = sum([
        abs(abs(minimized_fluxes[reaction_id]) - measured_flux)
        for reaction_id, measured_flux in measurements.items()
        if reaction_id in minimized_fluxes
    ])

    assert expected_distance == pytest.approx(solution.objective_value)
    assert expected_distance == pytest.approx(calculated_distance)
