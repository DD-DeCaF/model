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
from marshmallow.validate import OneOf

from simulations.modeling.gnomic_helpers import full_genotype


class Operation(Schema):
    operation = fields.String(required=True)
    type = fields.String(required=True)
    id = fields.String(missing=None)
    data = fields.Raw(missing=None)

    class Meta:
        strict = True


class MediumCompound(Schema):
    name = fields.String(required=True)
    # `id` should be a valid identifier in the namespace defined below.
    # `namespace` should match a namespace identifier from miriam.
    # See https://www.ebi.ac.uk/miriam/main/collections
    id = fields.String(required=True)
    namespace = fields.String(required=True)

    class Meta:
        strict = True


class Measurement(Schema):
    name = fields.String(required=True)
    # `id` should be a valid identifier in the namespace defined below.
    # `namespace` should match a namespace identifier from miriam.
    # See https://www.ebi.ac.uk/miriam/main/collections
    id = fields.String(required=True)
    namespace = fields.String(required=True)
    measurements = fields.List(fields.Float())
    type = fields.String(
        required=True, validate=OneOf(["compound", "reaction", "protein"])
    )

    class Meta:
        strict = True


class GrowthRate(Schema):
    measurements = fields.List(fields.Float())

    class Meta:
        strict = True


class ModificationRequest(Schema):
    medium = fields.Nested(MediumCompound, many=True, missing=[])
    genotype = fields.Function(deserialize=full_genotype, missing=None)
    growth_rate = fields.Nested(GrowthRate, missing=None)
    measurements = fields.Nested(Measurement, many=True, missing=[])

    class Meta:
        strict = True


class SimulationRequest(Schema):
    model_id = fields.Integer(required=True)
    method = fields.String(missing="fba")
    objective_id = fields.String(missing=None)
    objective_direction = fields.String(missing=None)
    operations = fields.Nested(Operation, many=True, missing=[])

    class Meta:
        strict = True
