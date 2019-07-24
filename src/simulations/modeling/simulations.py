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

import logging

from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis, pfba


logger = logging.getLogger(__name__)

METHODS = ["fba", "pfba", "fva", "pfba-fva"]


def simulate(model, biomass_reaction, method, objective_id, objective_direction):
    if method not in METHODS:
        raise ValueError(f"Unsupported simulation method '{method}'")

    if objective_id:
        model.objective = model.reactions.get_by_id(objective_id)
    if objective_direction:
        model.objective_direction = objective_direction

    try:
        logger.info(f"Simulating model {model.id} with {method}")
        if method == "fba":
            solution = model.optimize(raise_error=True)
        elif method == "pfba":
            solution = pfba(model)
        elif method == "fva":
            # FIXME: accept list of relevant fva reactions to calculate
            solution = flux_variability_analysis(model)
        elif method == "pfba-fva":
            # FIXME: accept list of relevant fva reactions to calculate
            solution = flux_variability_analysis(
                model, fraction_of_optimum=1, pfba_factor=1.05
            )
    except OptimizationError as error:
        logger.info(f"Optimization Error: {error}")
        flux_distribution = {}
        growth_rate = 0.0
    else:
        logger.info(f"Simulation completed successfully")
        if method in ("fba", "pfba"):
            flux_distribution = solution.fluxes.to_dict()
            growth_rate = flux_distribution[biomass_reaction]
        elif method in ("fva", "pfba-fva"):
            df = solution.rename(
                index=str, columns={"maximum": "upper_bound", "minimum": "lower_bound"}
            )
            for key in ["lower_bound", "upper_bound"]:
                df[key] = df[key].astype("float")
            flux_distribution = df.T.to_dict()
            growth_rate = flux_distribution[biomass_reaction]["upper_bound"]

    return flux_distribution, growth_rate
