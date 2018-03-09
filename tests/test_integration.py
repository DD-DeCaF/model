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
import pytest
from deepdiff import DeepDiff
import logging

from model.adapter import full_genotype, MediumChangeModel
from model.constants import METHODS, SIMULATION_METHOD, get_empty_changes
from model.storage import (restore_model, restore_from_db, save_changes_to_db, Models)
from model.operations import call_genes_to_reactions, modify_model, apply_reactions_add
from model.response import Response

logging.disable(logging.CRITICAL)

@pytest.mark.asyncio
async def test_call_genes_to_reactions():
    changes = full_genotype(['-aceA -sucCD +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'])
    result = await call_genes_to_reactions(changes)
    assert set(result.keys()) == {'BBa_J23100', 'AB326105', 'NP_600058', 'BBa_0010'}


@pytest.mark.asyncio
async def test_reactions_additions():
    ecoli_original = Models.get('iJO1366').copy()
    ecoli = ecoli_original.copy()
    ecoli.notes['changes'] = get_empty_changes()
    reactions = [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR81835', 'metabolites': None},
        {'id': 'MNXR83321', 'metabolites': None},
    ]
    reaction_ids = set([i['id'] for i in reactions])
    added_reactions = {'DM_12dgr182_9_12_e', 'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
                       'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
                       'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e',
                       'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    ecoli = await apply_reactions_add(ecoli, reactions)
    added_reactions_unique_ids = {i['id'] for i in ecoli.notes['changes']['added']['reactions']}
    assert len(ecoli.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    assert added_reactions_unique_ids - reaction_ids == added_reactions
    for reaction in ecoli.notes['changes']['added']['reactions']:
        assert ecoli.reactions.has_id(reaction['id'])
    reactions = [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR81835', 'metabolites': None},
    ]
    reaction_ids = set([i['id'] for i in reactions])
    ecoli = await apply_reactions_add(ecoli, reactions)
    added_reactions_unique_ids = {i['id'] for i in ecoli.notes['changes']['added']['reactions']}
    assert len(ecoli.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    assert added_reactions_unique_ids - reaction_ids == {
        'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
        'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e'
    }
    for reaction in ecoli.notes['changes']['added']['reactions']:
        assert ecoli.reactions.has_id(reaction['id'])
    removed_reactions = {'DM_12dgr182_9_12_e',
                         'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
                         'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    for reaction in removed_reactions:
        assert not ecoli.reactions.has_id(reaction)
    reactions = [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR81835', 'metabolites': None},
        {'id': 'MNXR83321', 'metabolites': None},
    ]
    reaction_ids = set([i['id'] for i in reactions])
    ecoli = await apply_reactions_add(ecoli, reactions)
    added_reactions_unique_ids = {i['id'] for i in ecoli.notes['changes']['added']['reactions']}
    assert len(ecoli.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    assert added_reactions_unique_ids - reaction_ids == added_reactions
    reaction_ids = {}
    ecoli = await apply_reactions_add(ecoli, list(reaction_ids))
    assert ecoli.notes['changes']['added']['reactions'] == []
    ecoli = await apply_reactions_add(ecoli, [{'id': 'MNXR83321', 'metabolites': None}, {'id': 'SUCR', 'metabolites': {'h2o_c': -1, 'sucr_c': -1, 'fru_c': 1, 'glc__D_c': 1}}])


@pytest.mark.asyncio
async def test_modify_model():
    message = {
        'to-return': ['tmy', 'fluxes', 'growth-rate', 'removed-reactions'],
        'objectives': ['chebi:17790'],
        'genotype-changes': ['+Aac'],
        'medium': [
            {'id': 'chebi:44080', 'concentration': 0.01},
            {'id': 'chebi:15075', 'concentration': 0.01},
            {'id': 'chebi:15377', 'concentration': 0.01},
            {'id': 'chebi:15378', 'concentration': 0.01},
            {'id': 'chebi:15379', 'concentration': 0.01},
            {'id': 'chebi:15982', 'concentration': 0.01},
            {'id': 'chebi:16189', 'concentration': 0.01},
            {'id': 'chebi:16526', 'concentration': 0.01},
            {'id': 'chebi:16643', 'concentration': 0.01},
            {'id': 'chebi:17883', 'concentration': 0.01},
            {'id': 'chebi:18212', 'concentration': 0.01},
            {'id': 'chebi:18367', 'concentration': 0.01},
            {'id': 'chebi:18420', 'concentration': 0.01},
            {'id': 'chebi:25371', 'concentration': 0.01},
            {'id': 'chebi:27638', 'concentration': 0.01},
            {'id': 'chebi:28938', 'concentration': 0.01},
            {'id': 'chebi:29033', 'concentration': 0.01},
            {'id': 'chebi:29034', 'concentration': 0.01},
            {'id': 'chebi:29035', 'concentration': 0.01},
            {'id': 'chebi:29036', 'concentration': 0.01},
            {'id': 'chebi:29101', 'concentration': 0.01},
            {'id': 'chebi:29103', 'concentration': 0.01},
            {'id': 'chebi:29105', 'concentration': 0.01},
            {'id': 'chebi:29108', 'concentration': 0.01},
            {'id': 'chebi:36271', 'concentration': 0.01},
            {'id': 'chebi:42758', 'concentration': 0.01},
            {'id': 'chebi:49786', 'concentration': 0.01}
        ],
        'measurements': [{'id': 'chebi:44080', 'measurements': [-15, -11, -14, -12], 'unit': 'mg', 'name': 'glucose',
                          'type': 'compound'},
                         {'id': 'PFK', 'measurements': [5, 5, 5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'}],
        'reactions-knockout': ['GLUDy', '3HAD160'],
    }
    wildtype = await restore_model('iJO1366')
    modified = await modify_model(message, wildtype.copy())
    assert 'EX_meoh_e' in modified.medium
    db_key = await save_changes_to_db(modified, 'iJO1366', message)
    restored_model = (await restore_from_db(db_key)).copy()
    assert restored_model.medium == modified.medium


FATTY_ACID_ECOLI = ['HACD2', 'ACACT1r', 'ECOAH3', 'HACD3', 'ECOAH1', 'ECOAH7', 'ACACT5r', 'ECOAH2', 'BUTCT',
                    'FACOAL60t2pp', 'ACACT6r', 'HACD7', 'ECOAH6', 'ACACT4r', 'FACOAL100t2pp', 'FACOAL160t2pp', 'ECOAH4',
                    'CTECOAI7', 'HACD1', 'ACOAD4f', 'FACOAL161t2pp', 'FACOAL80t2pp', 'CTECOAI8', 'HACD5', 'ACOAD1f',
                    'ACACT7r', 'ECOAH5', 'FACOAL180t2pp', 'CTECOAI6', 'FACOAL140t2pp', 'ACOAD7f', 'ACACT2r',
                    'FACOAL141t2pp', 'FACOAL181t2pp', 'HACD8', 'ACOAD8f', 'ACOAD2f', 'ECOAH8', 'ACACT3r', 'ACOAD3f',
                    'ACACT8r', 'HACD4', 'HACD6', 'ACOAD5f', 'ACOAD6f', 'FACOAL120t2pp']


@pytest.mark.asyncio
async def test_simulation_methods():
    for method in METHODS:
        message = {SIMULATION_METHOD: method}
        model = (await restore_model('iJO1366')).copy()
        response = Response(model, message)
        if method not in {'fva', 'pfba-fva'}:
            reactions_ids = [i.id for i in model.reactions]
            assert set(response.fluxes().keys()) == set(reactions_ids)


@pytest.mark.asyncio
async def test_restore_from_cache():
    wild_type_id = 'iMM904'
    message = {
        'to-return': ['model', 'fluxes'],
        'genotype-changes': [
            ('+promoter.Sc_TDH3:crtE_WT:terminator.Sc_CYC1,+promoter.Sc_TDH3:crtYB_WT:terminator.Sc_CYC1,'
             '+promoter.Sc_TDH3:crtI_WT:terminator.Sc_CYC1')],
        'simulation-method': 'pfba',
        'measurements': [{'unit': 'mmol', 'id': 'chebi:12965', 'measurements': [-1.1, -1.0],
                          'name': 'aldehydo-D-glucose', 'type': 'compound'},
                         {'unit': 'mmol', 'id': 'chebi:17579', 'measurements': [0.006, 0.008, 0.0065, 0.0007],
                          'name': 'beta-carotene', 'type': 'compound'},
                         {'id': 'PFK', 'measurements': [5], 'type': 'reaction', 'db_name': 'bigg.reaction'}]
    }
    model = await modify_model(message, (await restore_model(wild_type_id)).copy())
    db_key = await save_changes_to_db(model, wild_type_id, message)
    restored_model = (await restore_from_db(db_key)).copy()
    reactions = {
        r.id: dict(lower_bound=r.lower_bound, upper_bound=r.upper_bound)
        for r in model.reactions
    }
    reactions_restored = {
        r.id: dict(lower_bound=r.lower_bound, upper_bound=r.upper_bound)
        for r in restored_model.reactions
    }
    assert DeepDiff(reactions, reactions_restored) == {}
