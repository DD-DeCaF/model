import pytest
from deepdiff import DeepDiff
from model.app import call_genes_to_reactions, modify_model, restore_model, \
    METHODS, Response, SIMULATION_METHOD, FVA_REACTIONS, ProblemCache, \
    restore_from_db, save_changes_to_db
from driven.generic.adapter import full_genotype


@pytest.mark.asyncio
async def test_call_genes_to_reactions():
    changes = full_genotype(['-aceA -sucCD +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'])
    result = await call_genes_to_reactions(changes)
    assert set(result.keys()) == {'BBa_J23100', 'AB326105', 'NP_600058', 'BBa_0010'}


@pytest.mark.asyncio
async def test_modify_model():
    message = {
        'to-return': ['tmy', 'fluxes', 'growth-rate', 'removed-reactions'],
        'objectives': ['chebi:17790'],
        'genotype-changes': ['+Aac'],
        'medium': [{'id': 'chebi:44080', 'concentration': 0.01}],
        'measurements': [{'id': 'chebi:44080', 'measurement': -15, 'unit': 'mg', 'name': 'glucose'}],
        'reactions-knockout': ['GLUDy', '3HAD160'],
    }
    assert await modify_model(message, (await restore_model('iJO1366')).copy())


@pytest.mark.asyncio
async def test_simulation_methods():
    for method in METHODS:
        print(method)
        message = {SIMULATION_METHOD: method}
        if method == 'fva':
            message[FVA_REACTIONS] = ['MDH', 'ICL', 'nonsense']
        model = (await restore_model('iJO1366')).copy()
        cache = ProblemCache(model)
        response = Response(model, message, cache=cache)
        if method == 'fva':
            reactions_ids = ['MDH', 'ICL', 'BIOMASS_Ec_iJO1366_core_53p95M']
        else:
            reactions_ids = [i.id for i in model.reactions]
        assert set(response.fluxes().keys()) == set(reactions_ids)


@pytest.mark.asyncio
async def test_restore_from_cache():
    wild_type_id = 'iMM904'
    message = {
        'to-return': ['model', 'fluxes'],
        'genotype-changes': ['+promoter.Sc_TDH3:crtE_WT:terminator.Sc_CYC1,+promoter.Sc_TDH3:crtYB_WT:terminator.Sc_CYC1,+promoter.Sc_TDH3:crtI_WT:terminator.Sc_CYC1'],
        'simulation-method': 'pfba',
        'measurements': [{'unit': 'mmol', 'id': 'chebi:12965', 'measurement': -1.0, 'name': 'aldehydo-D-glucose'}, {'unit': 'mmol', 'id': 'chebi:17579', 'measurement': 0.0007, 'name': 'beta-carotene'}]
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
