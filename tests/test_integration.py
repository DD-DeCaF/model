import pytest
from model.app import call_genes_to_reactions, modify_model, restore_model
from driven.generic.adapter import full_genotype


@pytest.mark.asyncio
async def test_call_genes_to_reactions():
    changes = full_genotype(['-aceA -sucCD +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'])
    result = await call_genes_to_reactions(changes)
    assert set(result.keys()) == {'BBa_J23100', 'AB326105', 'NP_600058', 'BBa_0010'}


@pytest.mark.asyncio
async def test_modify_model():
    message = {
        'to-return': ['tmy', 'fluxes'],
        'objectives': ['chebi:17790'],
        'genotype-changes': ['+Aac'],
        'medium': [{'id': 'chebi:44080', 'concentration': 0.01}],
        'measurements': [{'id': 'chebi:44080', 'measurement': -15, 'unit': 'mg', 'name': 'glucose'}],
        'reactions-knockout': ['GLUDy', '3HAD160'],
    }
    assert await modify_model(message, restore_model('iJO1366'))
