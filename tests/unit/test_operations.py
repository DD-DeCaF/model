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

from model.operations import apply_operations


def test_add_reaction(e_coli_core):
    e_coli_core, biomass_reaction = e_coli_core
    assert not e_coli_core.reactions.has_id('FOOBAR')
    apply_operations(e_coli_core, [{
        'operation': "add",
        'type': "reaction",
        'data': {
            'id': 'FOOBAR',
            'name': 'Foo Bar',
            'metabolites': {
                'accoa_c': -1.0,
                'cit_c': 1.0,
                'coa_c': 1.0,
                'h2o_c': -1.0,
                'h_c': 1.0,
                'oaa_c': -1.0
            },
            'lower_bound': 0.0,
            'upper_bound': 1000.0,
            'gene_reaction_rule': 'b0720',
        }
    }])
    assert e_coli_core.reactions.FOOBAR.bounds == (0.0, 1000.0)


def test_modify_reaction(e_coli_core):
    e_coli_core, biomass_reaction = e_coli_core
    assert e_coli_core.reactions.CS.bounds == (0.0, 1000.0)
    apply_operations(e_coli_core, [{
        'operation': "modify",
        'type': "reaction",
        'id': "CS",
        'data': {
            'id': 'CS',
            'lower_bound': -20.0,
            'upper_bound': 20.0,
        }
    }])
    assert e_coli_core.reactions.CS.bounds == (-20.0, 20)


def test_knockout_reaction(e_coli_core):
    e_coli_core, biomass_reaction = e_coli_core
    assert e_coli_core.reactions.CS.bounds != (0.0, 0.0)
    apply_operations(e_coli_core, [{'operation': "knockout", 'type': "reaction", 'id': "CS"}])
    assert e_coli_core.reactions.CS.bounds == (0.0, 0.0)


def test_knockout_gene(e_coli_core):
    e_coli_core, biomass_reaction = e_coli_core
    assert e_coli_core.genes.b4025.functional
    assert all([r.bounds != (0.0, 0.0) for r in e_coli_core.genes.b4025.reactions])
    apply_operations(e_coli_core, [{'operation': "knockout", 'type': "gene", 'id': "b4025"}])
    assert not e_coli_core.genes.b4025.functional
    assert all([r.bounds == (0.0, 0.0) for r in e_coli_core.genes.b4025.reactions])
