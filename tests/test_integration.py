import pytest
from deepdiff import DeepDiff
from model.app import (call_genes_to_reactions, modify_model, restore_model,
                       METHODS, Response, SIMULATION_METHOD, ProblemCache,
                       restore_from_db, save_changes_to_db, find_in_memory, EMPTY_CHANGES, apply_reactions_add)
from model.adapter import full_genotype
from copy import deepcopy


@pytest.mark.asyncio
async def test_call_genes_to_reactions():
    changes = full_genotype(['-aceA -sucCD +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'])
    result = await call_genes_to_reactions(changes)
    assert set(result.keys()) == {'BBa_J23100', 'AB326105', 'NP_600058', 'BBa_0010'}


@pytest.mark.asyncio
async def test_reactions_additions():
    ecoli_original = find_in_memory('iJO1366').copy()
    ecoli = ecoli_original.copy()
    ecoli.notes['changes'] = deepcopy(EMPTY_CHANGES)
    reaction_ids = {'MNXR69355', 'MNXR81835', 'MNXR83321'}
    added_reactions = {'DM_12dgr182_9_12_e', 'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
            'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
            'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e',
            'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    ecoli = await apply_reactions_add(ecoli, list(reaction_ids))
    assert {i['id'] for i in ecoli.notes['changes']['added']['reactions']} - reaction_ids == added_reactions
    for reaction in ecoli.notes['changes']['added']['reactions']:
        assert ecoli.reactions.has_id(reaction['id'])
    reaction_ids -= {'MNXR83321'}
    ecoli = await apply_reactions_add(ecoli, list(reaction_ids))
    assert {i['id'] for i in ecoli.notes['changes']['added']['reactions']} - reaction_ids == \
           {'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
            'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e'}
    for reaction in ecoli.notes['changes']['added']['reactions']:
        assert ecoli.reactions.has_id(reaction['id'])
    removed_reactions = {'DM_12dgr182_9_12_e',
            'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
            'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    for reaction in removed_reactions:
        assert not ecoli.reactions.has_id(reaction)
    reaction_ids = {'MNXR69355', 'MNXR81835', 'MNXR83321'}
    ecoli = await apply_reactions_add(ecoli, list(reaction_ids))
    assert {i['id'] for i in ecoli.notes['changes']['added']['reactions']} - reaction_ids == added_reactions
    reaction_ids = {}
    ecoli = await apply_reactions_add(ecoli, list(reaction_ids))
    assert ecoli.notes['changes']['added']['reactions'] == []


@pytest.mark.asyncio
async def test_modify_model():
    message = {
        'to-return': ['tmy', 'fluxes', 'growth-rate', 'removed-reactions'],
        'objectives': ['chebi:17790'],
        'genotype-changes': ['+Aac'],
        'medium': [{'id': 'chebi:44080', 'concentration': 0.01}],
        'measurements': [{'id': 'chebi:44080', 'measurements': [-15, -11, -14, -12], 'unit': 'mg', 'name': 'glucose',
                          'type': 'compound'},
                         {'id': 'PFK', 'measurements': [5, 5, 5, 5], 'type': 'reaction'}],
        'reactions-knockout': ['GLUDy', '3HAD160'],
    }
    assert await modify_model(message, (await restore_model('iJO1366')).copy())


FATTY_ACID_ECOLI = ['HACD2', 'ACACT1r', 'ECOAH3', 'HACD3', 'ECOAH1', 'ECOAH7', 'ACACT5r', 'ECOAH2', 'BUTCT',
                    'FACOAL60t2pp', 'ACACT6r', 'HACD7', 'ECOAH6', 'ACACT4r', 'FACOAL100t2pp', 'FACOAL160t2pp', 'ECOAH4',
                    'CTECOAI7', 'HACD1', 'ACOAD4f', 'FACOAL161t2pp', 'FACOAL80t2pp', 'CTECOAI8', 'HACD5', 'ACOAD1f',
                    'ACACT7r', 'ECOAH5', 'FACOAL180t2pp', 'CTECOAI6', 'FACOAL140t2pp', 'ACOAD7f', 'ACACT2r',
                    'FACOAL141t2pp', 'FACOAL181t2pp', 'HACD8', 'ACOAD8f', 'ACOAD2f', 'ECOAH8', 'ACACT3r', 'ACOAD3f',
                    'ACACT8r', 'HACD4', 'HACD6', 'ACOAD5f', 'ACOAD6f', 'FACOAL120t2pp']


@pytest.mark.asyncio
async def test_simulation_methods():
    for method in METHODS:
        print(method)
        message = {SIMULATION_METHOD: method}
        model = (await restore_model('iJO1366')).copy()
        cache = ProblemCache(model)
        response = Response(model, message, cache=cache)
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
                         {'id': 'PFK', 'measurements': [5], 'type': 'reaction'}]
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
