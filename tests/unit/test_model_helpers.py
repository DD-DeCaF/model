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

from model.cobra_helpers import get_unique_metabolite
from model.exceptions import NoIDMapping


def test_existing_metabolite(iJO1366):
    assert get_unique_metabolite(iJO1366, 'chebi:17790') == get_unique_metabolite(
        iJO1366, 'meoh', db_name='bigg.metabolite')
    assert get_unique_metabolite(iJO1366, 'succ', db_name='bigg.metabolite').formula == 'C4H4O4'
    with pytest.raises(NoIDMapping):
        get_unique_metabolite(iJO1366, 'wrong_id')
