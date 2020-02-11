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

from simulations.modeling.reframed_helpers import generate_transactions


def test_generate_transactions_uptake():
    # Create mock exchanges
    id2name = {'M1': 'Metabolite 1'}
    exchanges_mock_more_uptake = {('A', 'M1'): 10, ('B', 'M1'): 5,
                                  ('C', 'M1'): -15, ('D', 'M1'): -5}

    result = generate_transactions(id2name, exchanges_mock_more_uptake)
    assert result[0][4]+result[1][4] == exchanges_mock_more_uptake[
        ('A', 'M1')];
    assert result[2][4] + result[3][4] == exchanges_mock_more_uptake[
        ('B', 'M1')];
    assert result[0][4]+result[2][4]+result[4][4] == - \
        exchanges_mock_more_uptake[('C', 'M1')];
    assert result[1][4]+result[3][4]+result[5][4] == - \
        exchanges_mock_more_uptake[('D', 'M1')];


def test_generate_transactions_secretion():
    # Create mock exchanges
    id2name = {'M1': 'Metabolite 1'}
    exchanges_mock_more_secretion = {('A', 'M1'): 20, ('B', 'M1'): 30,
                                     ('C', 'M1'): -15, ('D', 'M1'): -5}

    result = generate_transactions(id2name, exchanges_mock_more_secretion)
    assert result[0][4] + result[1][4] + result[4][4] == exchanges_mock_more_secretion[
        ('A', 'M1')];
    assert result[2[4] + result[3][4] + result[5][4] == \
           exchanges_mock_more_secretion[
               ('B', 'M1')];
    assert result[0][4] + result[2][4] == - \
    exchanges_mock_more_secretion[('C', 'M1')];
    assert result[1][4] + result[3][4] == - \
    exchanges_mock_more_secretion[('D', 'M1')];