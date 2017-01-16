import pytest
from model.app import call_genes_to_reactions, modify_model, restore_model, \
    METHODS, Response, SIMULATION_METHOD, FVA_REACTIONS, ProblemCache
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
            message[FVA_REACTIONS] = ['MDH', 'ICL']
        model = (await restore_model('iJO1366')).copy()
        cache = ProblemCache(model)
        response = Response(model, message, cache=cache)
        if method == 'fva':
            reactions_ids = ['MDH', 'ICL', 'BIOMASS_Ec_iJO1366_core_53p95M']
        else:
            reactions_ids = [i.id for i in model.reactions]
        assert set(response.fluxes().keys()) == set(reactions_ids)
