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


logger = logging.getLogger(__name__)


def apply_operations(model, operations):
    for operation in operations:
        if operation['operation'] == 'add' and operation['type'] == 'reaction':
            _add_reaction(model, operation['data'])
        elif operation['operation'] == 'modify' and operation['type'] == 'reaction':
            _modify_reaction(model, operation['id'], operation['data'])
        elif operation['operation'] == 'knockout' and operation['type'] == 'reaction':
            _knockout_reaction(model, operation['id'])
        elif operation['operation'] == 'knockout' and operation['type'] == 'gene':
            _knockout_gene(model, operation['id'])
        else:
            raise ValueError(f"Invalid operation: Cannot perform operation '{operation['operation']}' on type "
                             f"'{operation['type']}'")


def _parse_metabolite(metabolite_id):
    """Create a metabolite object from an identifier."""
    try:
        compartment = metabolite_id.rsplit("_", 1)[1]
    except IndexError:
        # Since we cannot parse the compartment from the identifier,
        # we assume it is in the cytosol.
        compartment = "c"
    return Metabolite(metabolite_id, compartment=compartment)


def _add_reaction(model, data):
    logger.debug(f"Adding reaction to model '{model.id}' from: {data}")
    metabolites = [
        _parse_metabolite(m) for m in data['metabolites']
        if m not in model.metabolites
    ]
    model.add_metabolites(metabolites)
    reaction = Reaction(
        id=data['id'],
        name=data['name'],
        lower_bound=data['lower_bound'],
        upper_bound=data['upper_bound']
    )
    model.add_reactions([reaction])
    reaction.add_metabolites(data['metabolites'])


def _modify_reaction(model, id, data):
    logger.debug(f"Modifying reaction '{id}' in model '{model.id}'")
    model.reactions.get_by_id(id).bounds = data['lower_bound'], data['upper_bound']


def _knockout_reaction(model, id):
    logger.debug(f"Removing reaction '{id}' from model '{model.id}'")
    model.reactions.get_by_id(id).knock_out()


def _knockout_gene(model, id):
    logger.debug(f"Removing gene '{id}' from model '{model.id}'")
    gene = model.genes.query(lambda g: id in (g.id, g.name))[0]
    gene.knock_out()
