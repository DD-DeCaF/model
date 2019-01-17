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

"""Marshmallow schemas for marshalling the API endpoints."""

from marshmallow import Schema, fields


class Operation(Schema):
    operation = fields.String(required=True)
    type = fields.String(required=True)
    id = fields.String()
    data = fields.Raw()

    class Meta:
        strict = True


class MediumCompound(Schema):
    id = fields.String(required=True)

    class Meta:
        strict = True


class Measurement(Schema):
    id = fields.String(required=True)
    measurements = fields.List(fields.Float())
    type = fields.String(required=True)

    class Meta:
        strict = True


class ModificationRequest(Schema):
    medium = fields.Nested(MediumCompound, many=True, missing=None)
    genotype = fields.List(fields.String(), missing=None)
    measurements = fields.Nested(Measurement, many=True, missing=None)

    class Meta:
        strict = True


class SimulationRequest(Schema):
    model_id = fields.String(missing=None)
    model = fields.Raw(missing=None)
    biomass_reaction = fields.String(missing=None)
    method = fields.String(missing='fba')
    objective_id = fields.String(missing=None)
    objective_direction = fields.String(missing=None)
    operations = fields.Nested(Operation, many=True, missing=[])

    class Meta:
        strict = True
