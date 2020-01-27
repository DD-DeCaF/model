# Copyright 2019 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

from os.path import dirname, join
# from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
from pytfa import ThermoModel
from pytfa.io import (
    load_thermoDB,
    apply_compartment_data,
    read_compartment_data,
)
from pytfa.optim import relax_dgo

from math import log

logger = logging.getLogger(__name__)


class HandlerThermo:
    """ Handler of pytfa.ThermoModel's.

    Exposes al the methods of regular cobra.Models when called, with the exception
    of `.tmfa()`, that fallbacks to the `.optimize()` method of pytfa.ThermoModel.

    The preparation and conversion to pytfa.ThermoModel are only performed if the
    tmfa() method is called, to avoid unnecesary overhead when metabolomics data
    is supplied but the simulations method is not TMFA.
    """

    def __init__(self, model, metabolomics=None):
        # TODO: check if the solver is properly handled by the package.
        self.cobra_model = model
        self.thermo_model = None
        # required as attribute since they are applied after
        # the conversion is produced (just called if .tmfa() is called)
        self.metabolomics = metabolomics

    def __getattr__(self, item):
        """Expose `cobra.Model` API."""
        return getattr(self.cobra_model, item)

    def _convert(self):
        """Extracts thermodynamic data from equilibrator to then create the
        Thermodynamic model."""
        # thermo_data = build_thermo_from_equilibrator(self.cobra_model)
        thermo_data = load_thermoDB(
            join(dirname(__file__), "../../../data/thermo_data.thermodb")
        )
        tmodel = ThermoModel(thermo_data, self.cobra_model)
        tmodel.solver = "optlang-glpk"
        # while not input from user, use some default data from iJ1366
        apply_compartment_data(
            tmodel,
            read_compartment_data(
                join(dirname(__file__), "../../../data/compartment_data.json")
            ),
        )
        # create the variables and the constraints for the TMFA LP problem
        tmodel.prepare()
        tmodel.convert()
        self.thermo_model = tmodel

    def apply_metabolomics(self):
        """Apply metabolomics measurements to the TMFA problem.

        It is important to note that uncertinty provided by the user is
        here applied and it is a different issue than thermodynamic uncertainty.
        """
        logger.debug(
            f"Aplying metabolomics data to ThermoModel '{self.cobra_model.id}'"
        )
        for measure in self.metabolomics:
            met = measure["identifier"]
            conc = measure["measurement"]
            the_conc_var = self.thermo_model.log_concentration.get_by_id(met)
            # Do not forget the variables in the model are logs !
            # uncertainty is very relevant in TMFA
            the_conc_var.ub = log(conc + measure["uncertainty"])
            the_conc_var.lb = log(conc - measure["uncertainty"])

    def tmfa(self):
        """Fallback to pytfa.ThermoModel `.optimize()` method.

        The optimization is done only after checking that the model has been converted.
        """
        logger.debug(f"Computing TMFA from ThermoModel {self.cobra_model}")
        if self.thermo_model is None:
            self._convert()
            if self.metabolomics:
                self.apply_metabolomics()
        return self.thermo_model.optimize()

    def safe_tmfa(self):
        """Relaxation of the THERMODYNAMICS (not metabolomics) constraints so
        the model can grow.

        TODO: account for input `growth_rate`. This should not be related to the
        relaxation in thermodynamics but in a further relaxation/minimization
        involving metabolomics measurementens (check the related built-in functions)
        """
        logger.warning(
            "Relaxation not implemented `.safe_tmfa()`, redirection"
            "to plain `.tmfa()``"
        )
        return self.tmfa()
        # NOTE: the problem with relax_dgo is that it performs tmfa innecesarily
        # several times, so it could be implemented with better perfomance here.
        solution = self.tmfa()
        if solution.value < self.tolerance:
            # this check should be done (and it is done) in the pytfa side;
            # however it does the relaxation anyways.
            self.thermo_model = relax_dgo(self.thermo_model.build_relaxation())[
                0
            ]
            solution = self.tmfa()
        return solution
