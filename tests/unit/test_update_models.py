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

import os

import pytest
from cobra.io import read_sbml_model

from model.modeling.cobra_helpers import find_metabolite
from tools.update_models import update_local_models


@pytest.mark.skip(reason="Local models are no longer stored in this service. Move logic and tests to model-warehouse")
def test_update_models(tmpdir):
    update_local_models('e_coli_core', tmpdir)
    model = read_sbml_model(os.path.join(tmpdir, 'e_coli_core.sbml.gz'))
    glucose = find_metabolite(model, 'CHEBI:42758', "CHEBI")
    assert glucose.id == 'glc__D_e'
    assert glucose.annotation['bigg.metabolite'] == 'glc__D'
