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

from math import isclose

from model.adapter import full_genotype
from model.constants import GENOTYPE_CHANGES, get_empty_changes
from model.operations import (
    add_reactions, apply_reactions_knockouts, build_string_from_metabolites, change_bounds,
    convert_measurements_to_mmol, convert_mg_to_mmol, modify_model, new_features_identifiers, phase_plane_to_dict,
    product_reaction_variable)


def is_close(x, y):
    return isclose(x, y, abs_tol=10e-6)


def test_mg_to_mmol():
    assert is_close(convert_mg_to_mmol(34.5, 58.4), 0.59075342)
    assert is_close(convert_mg_to_mmol(18, 18), 1.0)


def test_b_number(iJO1366):
    message = {
        GENOTYPE_CHANGES: ['-b3067,-b3172,-b1084'],
    }
    iJO1366 = modify_model(message, iJO1366)
    assert not iJO1366.genes.b3172.functional


def test_product_reaction_variable(e_coli_core):
    assert product_reaction_variable(e_coli_core, 'bigg:akg').id == 'EX_akg_e'
    assert product_reaction_variable(e_coli_core, 'bigg:e4p') is None


def test_phase_plane_to_dict(e_coli_core):
    result = phase_plane_to_dict(e_coli_core, 'bigg:glu__L')
    assert set(result.keys()) == {'EX_glu__L_e', 'objective_lower_bound', 'objective_upper_bound'}
    assert len(set([len(v) for v in result.values()])) == 1
    assert phase_plane_to_dict(e_coli_core, 'bigg:g3p') == {}


def test_new_features_identifiers():
    changes = full_genotype(['-A -B +promoter.C:#D:#E:terminator.F', '+G', '+Y -H'])
    result = new_features_identifiers(changes)
    assert set(result) == {'C', 'D', 'E', 'F', 'G', 'Y'}


def test_reactions_knockouts(iJO1366):
    iJO1366_copy = iJO1366.copy()
    iJO1366_copy.notes['changes'] = get_empty_changes()
    reaction_ids = {'GLUDy', 'GLUDy', '3HAD160'}
    GLUDy_upper_bound = iJO1366_copy.reactions.get_by_id('GLUDy').upper_bound
    assert GLUDy_upper_bound != 0
    iJO1366_copy = apply_reactions_knockouts(iJO1366_copy, list(reaction_ids))
    assert set([i['id'] for i in iJO1366_copy.notes['changes']['removed']['reactions']]) == reaction_ids
    assert iJO1366_copy.reactions.get_by_id('GLUDy').upper_bound == 0
    reaction_ids = reaction_ids - {'GLUDy'}
    iJO1366_copy = apply_reactions_knockouts(iJO1366_copy, list(reaction_ids))
    assert set([i['id'] for i in iJO1366_copy.notes['changes']['removed']['reactions']]) == {'3HAD160'}
    assert GLUDy_upper_bound == iJO1366_copy.reactions.get_by_id('GLUDy').upper_bound
    reaction_ids = reaction_ids - {'3HAD160'}
    iJO1366_copy = apply_reactions_knockouts(iJO1366_copy, list(reaction_ids))
    assert set([i['id'] for i in iJO1366_copy.notes['changes']['removed']['reactions']]) == set()
    assert is_close(iJO1366_copy.optimize().objective_value, iJO1366.optimize().objective_value)


def test_reactions_change_bounds(iJO1366):
    iJO1366_copy = iJO1366.copy()
    iJO1366_copy.notes['changes'] = get_empty_changes()
    reaction_ids = [{'id': "ACONTb", 'lower_bound': -3, 'upper_bound': 3.5},
                    {'id': "FBA3", 'lower_bound': -996, 'upper_bound': 1000}]
    FBA3_upper_bound = iJO1366_copy.reactions.get_by_id('FBA3').upper_bound
    assert FBA3_upper_bound != 0
    iJO1366_copy = change_bounds(iJO1366_copy, list(reaction_ids))
    reaction_ids = [{'id': "FBA3", 'lower_bound': -996, 'upper_bound': 1000}]
    iJO1366_copy = change_bounds(iJO1366_copy, list(reaction_ids))
    assert is_close(iJO1366_copy.optimize().objective_value, iJO1366.optimize().objective_value)


def test_convert_measurements_to_mmol(iJO1366):
    measurements = [{'id': 'chebi:17790', 'measurements': [32.04186], 'unit': 'mg', 'type': 'compound'}]
    assert convert_measurements_to_mmol(measurements, iJO1366) == [
        {'id': 'chebi:17790', 'measurements': [1.0], 'unit': 'mmol', 'type': 'compound'}]


def test_add_reactions(e_coli_core):
    keys = ('id', 'name', 'metabolites', 'lower_bound', 'upper_bound', 'gene_reaction_rule')
    info = ('newid', '', {}, 0, 0, '')
    changes = [dict(zip(keys, info))] * 2
    add_reactions(e_coli_core, changes)


def test_build_string_from_metabolites():
    metabolites = {'glc__D_c': -1, 'caro_c': 1}
    string = 'glc__D_c <=> caro_c'
    assert build_string_from_metabolites(metabolites) == string
