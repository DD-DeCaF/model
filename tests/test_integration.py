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

from deepdiff import DeepDiff

from model import adapter, storage
from model.adapter import full_genotype
from model.constants import get_empty_changes
from model.ice_client import ICE
from model.operations import apply_reactions_add, get_genotype_reactions, modify_model


def test_get_genotype_reactions():
    changes = full_genotype(['-aceA -sucCD +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'])
    result = get_genotype_reactions(changes)
    assert set(result.keys()) == {'BBa_J23100', 'AB326105', 'NP_600058', 'BBa_0010'}


def test_reactions_additions(monkeypatch, iJO1366):
    # Mock id-mapper api queries for efficiency
    def query_identifiers(object_ids, db_from, db_to):
        q = (object_ids, db_from, db_to)
        if q == (['MNXM2029', 'MNXM3447', 'MNXM368', 'MNXM7019'], 'mnx', 'chebi'):
            return {'MNXM368': ['16929', '10647', '12842', '26699', '52330', '57952']}
        elif q == (['MNXM2029', 'MNXM3447', 'MNXM368', 'MNXM7019'], 'mnx', 'bigg'):
            return {'MNXM3447': ['2agpe141'], 'MNXM368': ['g3pe'], 'MNXM7019': ['apg141'], 'MNXM2029': ['pg141']}
        elif q == (['MNXM1', 'MNXM2899', 'MNXM33', 'MNXM38', 'MNXM4923'], 'mnx', 'chebi'):
            return {'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM33': ['16238', '13315', '21125', '42388', '4956', '57692'], 'MNXM38': ['17877', '13316', '21126', '42427', '4957', '58307'], 'MNXM2899': ['10970', '22746', '3060', '57253'], 'MNXM4923': ['27639', '10948', '21124', '4731', '58519']}  # noqa
        elif q == (['MNXM1', 'MNXM2899', 'MNXM33', 'MNXM38', 'MNXM4923'], 'mnx', 'bigg'):
            return {'MNXM1': ['h'], 'MNXM33': ['fad'], 'MNXM38': ['fadh2'], 'MNXM2899': ['bzsuccoa'], 'MNXM4923': ['phitcoa']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
            return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
            return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
            return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
            return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
            return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
        elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
            return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
        elif q == (['h2o', 'sucr', 'fru', 'glc__D'], 'bigg', 'chebi'):
            return {'h2o': ['33813', '30490', '29412', '29375', '29356', '5594', '44641', '13419', '13365', '16234', '5585', '44819', '44701', '44292', '43228', '42857', '42043', '27313', '13352', '10743', '15377'], 'glc__D': ['17634', '12965', '20999', '4167'], 'fru': ['28757', '24104', '24110', '5172', '37721', '48095', '4119', '47424'], 'sucr': ['17992', '15128', '26812', '45795', '9314']}  # noqa
        elif q == (['h2o', 'sucr', 'fru', 'glc__D'], 'bigg', 'mnx'):
            return {'h2o': ['MNXM2'], 'glc__D': ['MNXM41'], 'fru': ['MNXM175'], 'sucr': ['MNXM167']}
        raise NotImplemented(f"Unmocked query!")
    monkeypatch.setattr(adapter, 'query_identifiers', query_identifiers)

    iJO1366.notes['changes'] = get_empty_changes()
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
    iJO1366 = apply_reactions_add(iJO1366, reactions)
    added_reactions_unique_ids = {i['id'] for i in iJO1366.notes['changes']['added']['reactions']}
    assert len(iJO1366.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    assert added_reactions_unique_ids - reaction_ids == added_reactions
    for reaction in iJO1366.notes['changes']['added']['reactions']:
        assert iJO1366.reactions.has_id(reaction['id'])
    reactions = [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR81835', 'metabolites': None},
    ]
    reaction_ids = set([i['id'] for i in reactions])
    iJO1366 = apply_reactions_add(iJO1366, reactions)
    added_reactions_unique_ids = {i['id'] for i in iJO1366.notes['changes']['added']['reactions']}
    assert len(iJO1366.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    assert added_reactions_unique_ids - reaction_ids == {
        'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
        'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e'
    }
    for reaction in iJO1366.notes['changes']['added']['reactions']:
        assert iJO1366.reactions.has_id(reaction['id'])
    removed_reactions = {'DM_12dgr182_9_12_e',
                         'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
                         'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    for reaction in removed_reactions:
        assert not iJO1366.reactions.has_id(reaction)
    reactions = [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR81835', 'metabolites': None},
        {'id': 'MNXR83321', 'metabolites': None},
    ]
    reaction_ids = set([i['id'] for i in reactions])
    iJO1366 = apply_reactions_add(iJO1366, reactions)
    added_reactions_unique_ids = {i['id'] for i in iJO1366.notes['changes']['added']['reactions']}
    assert len(iJO1366.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    assert added_reactions_unique_ids - reaction_ids == added_reactions
    reaction_ids = {}
    iJO1366 = apply_reactions_add(iJO1366, list(reaction_ids))
    assert iJO1366.notes['changes']['added']['reactions'] == []
    iJO1366 = apply_reactions_add(iJO1366, [{'id': 'MNXR83321', 'metabolites': None},
                                            {'id': 'SUCR', 'metabolites': {'h2o_c': -1,
                                                                           'sucr_c': -1,
                                                                           'fru_c': 1,
                                                                           'glc__D_c': 1}}])


def test_modify_model(iJO1366):
    message = {
        'to-return': ['tmy', 'fluxes', 'growth-rate', 'removed-reactions'],
        'theoretical-objectives': ['chebi:17790'],
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
    modified = modify_model(message, iJO1366)
    assert 'EX_meoh_e' in modified.medium
    db_key = storage.save_changes(modified, message)
    restored_model = storage.restore_from_key(db_key)
    assert restored_model.medium == modified.medium


FATTY_ACID_ECOLI = ['HACD2', 'ACACT1r', 'ECOAH3', 'HACD3', 'ECOAH1', 'ECOAH7', 'ACACT5r', 'ECOAH2', 'BUTCT',
                    'FACOAL60t2pp', 'ACACT6r', 'HACD7', 'ECOAH6', 'ACACT4r', 'FACOAL100t2pp', 'FACOAL160t2pp', 'ECOAH4',
                    'CTECOAI7', 'HACD1', 'ACOAD4f', 'FACOAL161t2pp', 'FACOAL80t2pp', 'CTECOAI8', 'HACD5', 'ACOAD1f',
                    'ACACT7r', 'ECOAH5', 'FACOAL180t2pp', 'CTECOAI6', 'FACOAL140t2pp', 'ACOAD7f', 'ACACT2r',
                    'FACOAL141t2pp', 'FACOAL181t2pp', 'HACD8', 'ACOAD8f', 'ACOAD2f', 'ECOAH8', 'ACACT3r', 'ACOAD3f',
                    'ACACT8r', 'HACD4', 'HACD6', 'ACOAD5f', 'ACOAD6f', 'FACOAL120t2pp']


def test_restore_from_cache(monkeypatch, iMM904):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, 'get_reaction_equations', lambda self, genotype: {})

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
    model = modify_model(message, iMM904)
    db_key = storage.save_changes(model, message)
    restored_model = storage.restore_from_key(db_key).copy()
    reactions = {
        r.id: dict(lower_bound=r.lower_bound, upper_bound=r.upper_bound)
        for r in model.reactions
    }
    reactions_restored = {
        r.id: dict(lower_bound=r.lower_bound, upper_bound=r.upper_bound)
        for r in restored_model.reactions
    }
    assert DeepDiff(reactions, reactions_restored) == {}
