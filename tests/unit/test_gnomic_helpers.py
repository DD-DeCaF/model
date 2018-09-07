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

from model.gnomic_helpers import full_genotype, feature_additions, feature_knockouts


def test_feature_operations():
    changes = full_genotype(['-A -B +promoter.C:#D:#E:terminator.F', '+G', '+B +Y -H'])
    assert set(feature_knockouts(changes)) == {'A', 'H'}
    assert set(feature_additions(changes)) == {'C', 'D', 'E', 'F', 'G', 'Y'}
