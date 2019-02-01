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

import pytest

from model.exceptions import MetaboliteNotFound
from model.modeling.cobra_helpers import find_metabolite


def test_existing_metabolite(iJO1366):
    iJO1366, biomass_reaction = iJO1366
    assert find_metabolite(iJO1366, 'CHEBI:17790', 'chebi', 'e') == find_metabolite(iJO1366, 'meoh', 'bigg.metabolite', 'e')
    assert find_metabolite(iJO1366, 'succ', 'bigg.metabolite', 'e').formula == 'C4H4O4'
    with pytest.raises(MetaboliteNotFound):
        find_metabolite(iJO1366, 'wrong_id', 'wrong_namespace', 'e')
