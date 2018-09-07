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

"""Test integrations between different modules"""

import pytest

from model import deltas
from model.adapter import adapt_from_medium
from model.operations import apply_operations


def test_modify_restore(iJO1366):
    """
    Modifying a model, saving it, loading the operations, and re-applying them to a wild type model, should yield a
    model with the same modifications.
    """
    original_medium = iJO1366.medium
    medium = [
        {'id': 'chebi:44080'},
        {'id': 'chebi:15075'},
        {'id': 'chebi:15377'},
        {'id': 'chebi:15378'},
    ]

    # Applying the new medium should yield a different composition
    with iJO1366:
        operations, errors = adapt_from_medium(iJO1366, medium)
        apply_operations(iJO1366, operations)
        modified_medium = iJO1366.medium
        assert original_medium != modified_medium

    # Cache and reload the operations
    deltas.save(iJO1366.id, {'medium': medium}, operations)
    loaded_operations = deltas.load(iJO1366.id, {'medium': medium})

    # Reapplying the loaded operations to a clean model should result in the exact same medium composition again
    with iJO1366:
        apply_operations(iJO1366, loaded_operations)
        assert iJO1366.medium == modified_medium


@pytest.mark.skip(reason="Namespaces are not currently mapped. Large test; unclear what exactly is tested for.")
def test_reactions_additions(monkeypatch, iJO1366):
    pass
    # # Mock id-mapper api queries for efficiency
    # def query_identifiers(object_ids, db_from, db_to):
    #     q = (object_ids, db_from, db_to)
    #     if q == (['MNXM2029', 'MNXM3447', 'MNXM368', 'MNXM7019'], 'mnx', 'chebi'):
    #         return {'MNXM368': ['16929', '10647', '12842', '26699', '52330', '57952']}
    #     elif q == (['MNXM2029', 'MNXM3447', 'MNXM368', 'MNXM7019'], 'mnx', 'bigg'):
    #         return {'MNXM3447': ['2agpe141'], 'MNXM368': ['g3pe'], 'MNXM7019': ['apg141'], 'MNXM2029': ['pg141']}
    #     elif q == (['MNXM1', 'MNXM2899', 'MNXM33', 'MNXM38', 'MNXM4923'], 'mnx', 'chebi'):
    #         return {'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM33': ['16238', '13315', '21125', '42388', '4956', '57692'], 'MNXM38': ['17877', '13316', '21126', '42427', '4957', '58307'], 'MNXM2899': ['10970', '22746', '3060', '57253'], 'MNXM4923': ['27639', '10948', '21124', '4731', '58519']}  # noqa
    #     elif q == (['MNXM1', 'MNXM2899', 'MNXM33', 'MNXM38', 'MNXM4923'], 'mnx', 'bigg'):
    #         return {'MNXM1': ['h'], 'MNXM33': ['fad'], 'MNXM38': ['fadh2'], 'MNXM2899': ['bzsuccoa'], 'MNXM4923': ['phitcoa']}  # noqa
    #     elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
    #         return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
    #     elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
    #         return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
    #     elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
    #         return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
    #     elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
    #         return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
    #     elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'chebi'):
    #         return {'MNXM89795': ['18307', '13487', '13495', '22100', '42751', '9811', '58439', '66914', '67119'], 'MNXM1': ['24636', '5584', '13357', '10744', '15378'], 'MNXM17': ['17659', '13445', '27230', '46402', '9802', '58223']}  # noqa
    #     elif q == (['MNXM147347', 'MNXM89795', 'MNXM1', 'MNXM146474', 'MNXM17'], 'mnx', 'bigg'):
    #         return {'MNXM147347': ['12dgr182_9_12'], 'MNXM146474': ['mgdg182_9_12'], 'MNXM89795': ['udpgal'], 'MNXM1': ['h'], 'MNXM17': ['udp']}  # noqa
    #     elif q == (['h2o', 'sucr', 'fru', 'glc__D'], 'bigg', 'chebi'):
    #         return {'h2o': ['33813', '30490', '29412', '29375', '29356', '5594', '44641', '13419', '13365', '16234', '5585', '44819', '44701', '44292', '43228', '42857', '42043', '27313', '13352', '10743', '15377'], 'glc__D': ['17634', '12965', '20999', '4167'], 'fru': ['28757', '24104', '24110', '5172', '37721', '48095', '4119', '47424'], 'sucr': ['17992', '15128', '26812', '45795', '9314']}  # noqa
    #     elif q == (['h2o', 'sucr', 'fru', 'glc__D'], 'bigg', 'mnx'):
    #         return {'h2o': ['MNXM2'], 'glc__D': ['MNXM41'], 'fru': ['MNXM175'], 'sucr': ['MNXM167']}
    #     raise NotImplemented(f"Unmocked query!")
    # monkeypatch.setattr(adapter, 'query_identifiers', query_identifiers)

    # iJO1366.notes['changes'] = get_empty_changes()
    # reactions = [
    #     {'id': 'MNXR69355', 'metabolites': None},
    #     {'id': 'MNXR81835', 'metabolites': None},
    #     {'id': 'MNXR83321', 'metabolites': None},
    # ]
    # reaction_ids = set([i['id'] for i in reactions])
    # added_reactions = {'DM_12dgr182_9_12_e', 'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
    #                    'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
    #                    'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e',
    #                    'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    # iJO1366 = apply_reactions_add(iJO1366, reactions)
    # added_reactions_unique_ids = {i['id'] for i in iJO1366.notes['changes']['added']['reactions']}
    # assert len(iJO1366.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    # assert added_reactions_unique_ids - reaction_ids == added_reactions
    # for reaction in iJO1366.notes['changes']['added']['reactions']:
    #     assert iJO1366.reactions.has_id(reaction['id'])
    # reactions = [
    #     {'id': 'MNXR69355', 'metabolites': None},
    #     {'id': 'MNXR81835', 'metabolites': None},
    # ]
    # reaction_ids = set([i['id'] for i in reactions])
    # iJO1366 = apply_reactions_add(iJO1366, reactions)
    # added_reactions_unique_ids = {i['id'] for i in iJO1366.notes['changes']['added']['reactions']}
    # assert len(iJO1366.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    # assert added_reactions_unique_ids - reaction_ids == {
    #     'DM_phitcoa_e', 'adapter_bzsuccoa_c_bzsuccoa_e',
    #     'adapter_phitcoa_c_phitcoa_e', 'DM_bzsuccoa_e'
    # }
    # for reaction in iJO1366.notes['changes']['added']['reactions']:
    #     assert iJO1366.reactions.has_id(reaction['id'])
    # removed_reactions = {'DM_12dgr182_9_12_e',
    #                      'DM_mgdg182_9_12_e', 'adapter_mgdg182_9_12_c_mgdg182_9_12_e',
    #                      'adapter_12dgr182_9_12_c_12dgr182_9_12_e'}
    # for reaction in removed_reactions:
    #     assert not iJO1366.reactions.has_id(reaction)
    # reactions = [
    #     {'id': 'MNXR69355', 'metabolites': None},
    #     {'id': 'MNXR81835', 'metabolites': None},
    #     {'id': 'MNXR83321', 'metabolites': None},
    # ]
    # reaction_ids = set([i['id'] for i in reactions])
    # iJO1366 = apply_reactions_add(iJO1366, reactions)
    # added_reactions_unique_ids = {i['id'] for i in iJO1366.notes['changes']['added']['reactions']}
    # assert len(iJO1366.notes['changes']['added']['reactions']) == len(added_reactions_unique_ids)
    # assert added_reactions_unique_ids - reaction_ids == added_reactions
    # reaction_ids = {}
    # iJO1366 = apply_reactions_add(iJO1366, list(reaction_ids))
    # assert iJO1366.notes['changes']['added']['reactions'] == []
    # iJO1366 = apply_reactions_add(iJO1366, [{'id': 'MNXR83321', 'metabolites': None},
    #                                         {'id': 'SUCR', 'metabolites': {'h2o_c': -1,
    #                                                                        'sucr_c': -1,
    #                                                                        'fru_c': 1,
    #                                                                        'glc__D_c': 1}}])
