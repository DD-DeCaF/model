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
import requests
from cobra import Model
from flask import g

from simulations import storage
from simulations.exceptions import Forbidden, Unauthorized


class MockResponseSuccess:
    status_code = 200

    def json(self):
        return {
            "id": 1,
            "model_serialized": {
                "version": "1",
                "id": "foo",
                "genes": [],
                "reactions": [],
                "metabolites": [],
                "compartments": {},
            },
            "project_id": None,
            "organism_id": 2,
            "default_biomass_reaction": "baz",
            "ec_model": False,
        }

    def raise_for_status(self):
        pass


class MockResponseUnauthorized:
    status_code = 401

    def json(self):
        return {"message": "Mocked response: Unauthorized"}


class MockResponseForbidden:
    status_code = 403

    def json(self):
        return {"message": "Mocked response: Forbidden"}


def test_get_model(monkeypatch, app):
    monkeypatch.setattr(requests, "get", lambda url, headers: MockResponseSuccess())
    g.jwt_valid = False
    assert type(storage.get(10).model) == Model


def test_get_model_forbidden(monkeypatch, app):
    monkeypatch.setattr(requests, "get", lambda url, headers: MockResponseForbidden())
    g.jwt_valid = False
    with pytest.raises(Forbidden):
        storage.get(11)


def test_get_model_unauthorized(monkeypatch, app):
    monkeypatch.setattr(
        requests, "get", lambda url, headers: MockResponseUnauthorized()
    )
    g.jwt_valid = False
    with pytest.raises(Unauthorized):
        storage.get(11)
