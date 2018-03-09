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
from math import isclose
import pytest

from model.adapter import full_genotype
from model.constants import GENOTYPE_CHANGES, get_empty_changes
from model.operations import (convert_mg_to_mmol, modify_model, product_reaction_variable,
                              phase_plane_to_dict, new_features_identifiers, add_reactions,
                              apply_reactions_knockouts, convert_measurements_to_mmol,
                              build_string_from_metabolites)
from model.storage import restore_model, Models

logging.disable(logging.CRITICAL)


def is_close(x, y):
    return isclose(x, y, abs_tol=10e-6)


def test_mg_to_mmol():
    assert is_close(convert_mg_to_mmol(34.5, 58.4), 0.59075342)
    assert is_close(convert_mg_to_mmol(18, 18), 1.0)


@pytest.mark.asyncio
async def test_b_number():
    model_id = 'iJO1366'
    message = {
        GENOTYPE_CHANGES: ['-b3067,-b3172,-b1084'],
    }
    model = await modify_model(message, (await restore_model(model_id)).copy())
    assert not model.genes.b3172.functional


def test_product_reaction_variable():
    ecoli = Models.get('iJO1366')
    assert product_reaction_variable(ecoli, 'bigg:akg').id == 'EX_akg_e'
    assert product_reaction_variable(ecoli, 'bigg:e4p') is None


def test_phase_plane_to_dict():
    ecoli = Models.get('iJO1366')
    result = phase_plane_to_dict(ecoli, 'bigg:glu__L')
    assert set(result.keys()) == {'EX_glu__L_e', 'objective_lower_bound', 'objective_upper_bound'}
    assert len(set([len(v) for v in result.values()])) == 1
    assert phase_plane_to_dict(ecoli, 'bigg:g3p') == {}


def test_new_features_identifiers():
    changes = full_genotype(['-A -B +promoter.C:#D:#E:terminator.F', '+G', '+Y -H'])
    result = new_features_identifiers(changes)
    assert set(result) == {'C', 'D', 'E', 'F', 'G', 'Y'}


@pytest.mark.asyncio
async def test_reactions_knockouts():
    ecoli_original = Models.get('iJO1366').copy()
    ecoli = ecoli_original.copy()
    ecoli.notes['changes'] = get_empty_changes()
    reaction_ids = {'GLUDy', 'GLUDy', '3HAD160'}
    GLUDy_upper_bound = ecoli.reactions.get_by_id('GLUDy').upper_bound
    assert GLUDy_upper_bound != 0
    ecoli = await apply_reactions_knockouts(ecoli, list(reaction_ids))
    assert set([i['id'] for i in ecoli.notes['changes']['removed']['reactions']]) == reaction_ids
    assert ecoli.reactions.get_by_id('GLUDy').upper_bound == 0
    reaction_ids = reaction_ids - {'GLUDy'}
    ecoli = await apply_reactions_knockouts(ecoli, list(reaction_ids))
    assert set([i['id'] for i in ecoli.notes['changes']['removed']['reactions']]) == {'3HAD160'}
    assert GLUDy_upper_bound == ecoli.reactions.get_by_id('GLUDy').upper_bound
    reaction_ids = reaction_ids - {'3HAD160'}
    ecoli = await apply_reactions_knockouts(ecoli, list(reaction_ids))
    assert set([i['id'] for i in ecoli.notes['changes']['removed']['reactions']]) == set()
    assert is_close(ecoli.optimize().objective_value, ecoli_original.optimize().objective_value)


@pytest.mark.asyncio
async def test_convert_measurements_to_mmol():
    ecoli = Models.get('iJO1366').copy()
    measurements = [{'id': 'chebi:17790', 'measurements': [32.04186], 'unit': 'mg', 'type': 'compound'}]
    assert convert_measurements_to_mmol(measurements, ecoli) == [
        {'id': 'chebi:17790', 'measurements': [1.0], 'unit': 'mmol', 'type': 'compound'}]


def test_add_reactions():
    ecoli = Models.get('iJO1366').copy()
    keys = ('id', 'name', 'metabolites', 'lower_bound', 'upper_bound', 'gene_reaction_rule')
    info = ('newid', '', {}, 0, 0, '')
    changes = [dict(zip(keys, info))] * 2
    add_reactions(ecoli, changes)


def test_build_string_from_metabolites():
    metabolites = {'glc__D_c': -1, 'caro_c': 1}
    string = 'glc__D_c <=> caro_c'
    assert build_string_from_metabolites(metabolites) == string
