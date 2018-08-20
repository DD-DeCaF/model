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

import random

from model import deltas


def test_delta_key_from_conditions():
    model_id = 'model_id'
    for _ in range(50):
        values1 = random.sample(range(20), 10)
        values2 = random.sample(range(20), 10)
        values3 = random.shuffle(values2[:])
        conditions1 = {'measurements': values1, 'genotype': values2}
        conditions2 = {'measurements': values2}
        conditions3 = {'measurements': values3, 'genotype': values2}
        conditions4 = {'genotype': values2, 'measurements': values1}
        assert deltas._delta_key(model_id, conditions1) != deltas._delta_key(model_id, conditions2)
        assert deltas._delta_key(model_id, conditions1) != deltas._delta_key(model_id, conditions3)
        assert deltas._delta_key(model_id, conditions1) == deltas._delta_key(model_id, conditions4)


def test_save_and_load(client):
    operations = [{'operation': 'remove', 'type': 'reaction', 'id': 'EX_ca2_e'}]
    conditions = {'medium': [{'id': 'chebi:42758'}]}
    key = deltas.save('model_foo', conditions, operations)
    assert deltas.load('model_foo', conditions) == operations
    assert deltas.load_from_key(key) == operations
