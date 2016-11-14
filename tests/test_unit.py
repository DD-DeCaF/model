import pytest
from model.app import existing_metabolite, NoIDMapping, restore_model, product_reaction_variable, phase_plane_to_dict, \
    new_features_identifiers, apply_reactions_knockouts, respond
from driven.generic.adapter import full_genotype


def test_existing_metabolite():
    ecoli = restore_model('iJO1366')
    assert existing_metabolite(ecoli, 'chebi:17790') == existing_metabolite(ecoli, 'bigg:meoh')
    assert existing_metabolite(ecoli, 'bigg:succ').formula == 'C4H4O4'
    with pytest.raises(NoIDMapping):
        existing_metabolite(ecoli, 'wrong_id')


def test_product_reaction_variable():
    ecoli = restore_model('iJO1366')
    assert product_reaction_variable(ecoli, 'bigg:akg').id == 'EX_akg_e'
    assert product_reaction_variable(ecoli, 'bigg:e4p') is None


def test_phase_plane_to_dict():
    ecoli = restore_model('iJO1366')
    result = phase_plane_to_dict(ecoli, 'bigg:glu__L')
    assert set(result.keys()) == {'EX_glu__L_e', 'objective_lower_bound', 'objective_upper_bound'}
    assert len(set([len(v) for v in result.values()])) == 1
    assert phase_plane_to_dict(ecoli, 'bigg:g3p') == {}


def test_new_features_identifiers():
    changes = full_genotype(['-A -B +promoter.C:#D:#E:terminator.F', '+G', '+Y -H'])
    result = new_features_identifiers(changes)
    assert set(result) == {'C', 'D', 'E', 'F', 'G', 'Y'}


def test_respond():
    message = {'to-return': ['fluxes', 'tmy', 'model'], 'objectives': ['bigg:akg']}
    assert set(respond(message, restore_model('iJO1366')).keys()) == set(message['to-return'])


@pytest.mark.asyncio
async def test_apply_reactions_knockouts():
    ecoli = restore_model('iJO1366')
    result = await apply_reactions_knockouts(ecoli, ['GLUDy', '3HAD160', 'GLUDy'])
    assert result.reactions.GLUDy.lower_bound == result.reactions.GLUDy.upper_bound == 0
