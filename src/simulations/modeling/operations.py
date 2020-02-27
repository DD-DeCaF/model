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

from cobra import Metabolite, Reaction

from simulations.exceptions import CompartmentNotFound
from simulations.modeling.cobra_helpers import parse_bigg_compartment
from simulations.modeling.pytfa_helpers import HandlerThermo


logger = logging.getLogger(__name__)


def apply_operations(model, operations):
    for operation in operations:
        if operation["operation"] == "add" and operation["type"] == "reaction":
            _add_reaction(model, operation["data"])
        elif operation["operation"] == "modify" and operation["type"] == "reaction":
            _modify_reaction(model, operation["id"], operation["data"])
        elif operation["operation"] == "modify" and operation["type"] == "thermo":
            _apply_metabolomics(model, operation["data"])
        elif operation["operation"] == "knockout" and operation["type"] == "reaction":
            _knockout_reaction(model, operation["id"])
        elif operation["operation"] == "knockout" and operation["type"] == "gene":
            _knockout_gene(model, operation["id"])
        elif operation["operation"] == "remove" and operation["type"] == "reaction":
            _remove_reaction(model, operation["id"])
        else:
            raise ValueError(
                f"Invalid operation: Cannot perform operation "
                f"'{operation['operation']}' on type '{operation['type']}'"
            )


def _parse_metabolite(metabolite_id, model):
    """Create a metabolite object from an identifier."""
    try:
        _, compartment_id = parse_bigg_compartment(metabolite_id, model)
    except (ValueError, CompartmentNotFound):
        # Since we cannot parse the compartment from the identifier,
        # we assume it is in the cytosol.
        compartment_id = "c"
    return Metabolite(metabolite_id, compartment=compartment_id)


def _add_reaction(model, data):
    logger.debug(f"Adding reaction to model '{model.id}' from: {data}")
    metabolites = [
        _parse_metabolite(m, model)
        for m in data["metabolites"]
        if m not in model.metabolites
    ]
    model.add_metabolites(metabolites)
    reaction = Reaction(
        id=data["id"],
        name=data["name"],
        lower_bound=data["lower_bound"],
        upper_bound=data["upper_bound"],
    )
    model.add_reactions([reaction])
    reaction.add_metabolites(data["metabolites"])


def _modify_reaction(model, id, data):
    logger.debug(
        f"Setting bounds of reaction '{id}' to ({data['lower_bound']}, "
        f"{data['upper_bound']}) in model '{model.id}'"
    )
    model.reactions.get_by_id(id).bounds = data["lower_bound"], data["upper_bound"]


def _knockout_reaction(model, id):
    logger.debug(f"Knocking out reaction '{id}' in model '{model.id}'")
    model.reactions.get_by_id(id).knock_out()


def _remove_reaction(model, id):
    logger.debug(f"Removing reaction '{id}' from model '{model.id}'")
    model.remove_reactions([model.reactions.get_by_id(id)])


def _knockout_gene(model, id):
    logger.debug(f"Knocking out gene '{id}' in model '{model.id}'")
    gene = model.genes.query(lambda g: id in (g.id, g.name))[0]
    gene.knock_out()


def _apply_metabolomics(model, metabolomics):
    logger.debug(f"Creating thermodynamic model from '{model.id}'")
    model = HandlerThermo(model, metabolomics)
