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
from cobra.io import read_sbml_model

from simulations import storage
from simulations.app import app as app_
from simulations.app import init_app


@pytest.fixture(scope="session")
def app():
    """Provide the initialized Flask app"""
    init_app(app_, None)
    with app_.app_context():
        yield app_


@pytest.fixture(scope="session")
def client(app):
    """Provide a Flask test client"""
    with app.test_client() as client:
        yield client


@pytest.fixture(scope="session")
def models():
    """
    Preload the storage module with test models.

    This fixture ensures models are loaded locally, and only once per test session.
    The returned dict contains a map of model identifier to the corresponding numeric id
    in the storage module.

    For endpoint tests, use the returned dict to access the numeric id for a given model
    id. For unit tests, consider using the function-scoped fixtures below to be able to
    make revertable modifications to the models.
    """
    model_keys = {
        "e_coli_core": 1,
        "e_coli_core_proprietary": 2,
        "iJO1366": 3,
        "eciML1515": 4,
    }

    model = read_sbml_model("tests/data/e_coli_core.sbml.gz")
    storage._MODELS[model_keys["e_coli_core"]] = storage.ModelWrapper(
        model, None, "Escherichia coli", "BIOMASS_Ecoli_core_w_GAM", False
    )
    model = read_sbml_model("tests/data/e_coli_core.sbml.gz")
    storage._MODELS[model_keys["e_coli_core_proprietary"]] = storage.ModelWrapper(
        model, 1, "Escherichia coli", "BIOMASS_Ecoli_core_w_GAM", False
    )
    model = read_sbml_model("tests/data/iJO1366.sbml.gz")
    storage._MODELS[model_keys["iJO1366"]] = storage.ModelWrapper(
        model, None, "Escherichia coli", "BIOMASS_Ec_iJO1366_core_53p95M", False
    )
    model = read_sbml_model("tests/data/eciML1515.xml.gz")
    storage._MODELS[model_keys["eciML1515"]] = storage.ModelWrapper(
        model, None, "Escherichia coli", "BIOMASS_Ec_iML1515_core_75p37M", True
    )
    return model_keys


@pytest.fixture(scope="function")
def e_coli_core(models):
    """
    Provide the e_coli_core model in a context manager, so that modifications are not
    persisted beyond the scope of the test function. This model is fairly small and
    should be preferred in test cases where possible.
    """
    wrapper = storage._MODELS[models["e_coli_core"]]
    with wrapper.model as model:
        yield model, wrapper.biomass_reaction, wrapper.is_ec_model


@pytest.fixture(scope="function")
def iJO1366(models):
    """
    Provide the iJO1366 model in a context manager, so that modifications are not
    persisted beyond the scope of the test function.
    """
    wrapper = storage._MODELS[models["iJO1366"]]
    with wrapper.model as model:
        yield model, wrapper.biomass_reaction, wrapper.is_ec_model


@pytest.fixture(scope="function")
def eciML1515(models):
    """
    Provide the eciML1515 model in a context manager, so that modifications are not
    persisted beyond the scope of the test function. This is an enzyme-constrained
    model, so it should only be used for testing enzyme-constrained related functions
    (e.g. direct integration of proteomics data).
    """
    wrapper = storage._MODELS[models["eciML1515"]]
    with wrapper.model as model:
        yield model, wrapper.biomass_reaction, wrapper.is_ec_model
