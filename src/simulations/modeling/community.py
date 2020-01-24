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
import tempfile

import cobra
import reframed


logger = logging.getLogger(__name__)


def simulate(models, medium):
    """
    Run a SteadyCom community simulation.

    Parameters
    ----------
    models: list(cobra.Model)
        A list of cobrapy model instances.
    medium: list(str)
        A list of compound names. Exchange reaction identifiers are assumed to
        be formatted according to: "EX_{compound}_e"
    """
    logger.debug("Converting cobrapy models to reframed models")
    rf_models = []
    for model in models:
        # The most funcational approach (albeit slow) seems to be to write and
        # reload SBML. reframed's cobrapy integration is currently pretty
        # minimal.
        with tempfile.NamedTemporaryFile() as file_:
            cobra.io.write_sbml_model(model, file_.name)
            rf_models.append(reframed.load_cbmodel(file_.name))

    logger.debug("Merging individual models to a community")
    community = reframed.Community("community", rf_models)

    logger.debug("Applying medium to the community")
    environment = reframed.Environment.from_compounds(
        medium, fmt_func=lambda x: f"R_EX_M_{x}_e"
    )
    environment.apply(community.merged_model, inplace=True)

    logger.info(f"Simulating community model with SteadyCom")
    solution = reframed.SteadyCom(community)

    logger.debug(f"Formatting solution response")
    return {
        "growth_rate": solution.growth,
        "abundance": solution.abundance,
        "cross_feeding": solution.cross_feeding(),
    }
