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

from cobra.io.dict import reaction_from_dict


logger = logging.getLogger(__name__)


def apply_operations(model, operations):
    for operation in operations:
        if operation['operation'] == 'add':
            if operation['type'] == 'reaction':
                _add_reaction(model, operation['id'], operation['data'])
            else:
                raise ValueError(f"Invalid operation: Cannot add type '{operation['type']}'")
        elif operation['operation'] == 'modify':
            if operation['type'] == 'reaction':
                _modify_reaction(model, operation['id'], operation['data'])
            else:
                raise ValueError(f"Invalid operation: Cannot modify type '{operation['type']}'")
        elif operation['operation'] == 'remove':
            if operation['type'] == 'reaction':
                _remove_reaction(model, operation['id'])
            elif operation['type'] == 'gene':
                _remove_gene(model, operation['id'])
            else:
                raise ValueError(f"Invalid operation: Cannot modify type '{operation['type']}'")
        else:
            raise ValueError(f"Invalid operation: Cannot perform operation '{operation['operation']}'")


def _add_reaction(model, id, data):
    logger.debug(f"Adding reaction '{id}' to model '{model.id}'")
    reaction = reaction_from_dict(data)
    model.add_reaction(reaction)


def _modify_reaction(model, id, data):
    logger.debug(f"Modifying reaction '{id}' in model '{model.id}'")
    model.reactions.get_by_id(id).bounds = data['bounds']


def _remove_reaction(model, id):
    logger.debug(f"Removing reaction '{id}' from model '{model.id}'")
    model.reactions.get_by_id(id).knock_out()


def _remove_gene(model, id):
    logger.debug(f"Removing gene '{id}' from model '{model.id}'")
    # NOTES(Ali): always looking for name+id in genes now.
    # NOTES(Ali): previous code: gene = model.genes.query(gene['name'], attribute="name")[0]
    # NOTES(Ali): clarify whether that is necessary
    gene = model.genes.query(lambda g: id in (g.id, g.name))[0]
    gene.knock_out()
