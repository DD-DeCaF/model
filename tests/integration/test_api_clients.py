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

"""Test external API service integrations"""

import pytest

from model.exceptions import PartNotFound
from model.ice_client import ICE


ice = ICE()


@pytest.mark.skip(reason="ICE seems to be occasionally unresponsive and halts CI builds")
@pytest.mark.parametrize("part", ('NP_600058', 'BBa_J23100', 'AB326105'))
def test_ice_missing_parts(part):
    """
    Query ICE for missing parts.
    Prone to failure when live data changes and will need to be kept in sync accordingly.
    """
    with pytest.raises(PartNotFound):
        ice.get_reaction_equations(part)


@pytest.mark.skip(reason="ICE seems to be occasionally unresponsive and halts CI builds")
def test_ice_existing_part():
    """
    Query ICE for an existing part.
    Prone to failure when live data changes and will need to be kept in sync accordingly.
    """
    result = ice.get_reaction_equations('BBa_0010')
    assert result == {'DECARB': 'acon_C <=> itacon + co2'}
