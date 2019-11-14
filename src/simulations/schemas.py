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

import gnomic
from marshmallow import Schema, fields


# For all reaction and compound references: `namespace` should match a namespace
# identifier from miriam[1] and `identifier` should be a valid identifier in that
# namespace.
# [1] https://www.ebi.ac.uk/miriam/main/collections


class Operation(Schema):
    operation = fields.String(required=True)
    type = fields.String(required=True)
    id = fields.String(missing=None)
    data = fields.Raw(missing=None)


class MediumCompound(Schema):
    name = fields.String(required=True)
    identifier = fields.String(required=True)
    namespace = fields.String(required=True)
    mass_concentration = fields.Float(required=True, allow_none=True)


class Fluxomics(Schema):
    name = fields.String(required=True)
    identifier = fields.String(required=True)
    namespace = fields.String(required=True)
    measurement = fields.Float(required=True)
    uncertainty = fields.Float(required=True)


class Metabolomics(Schema):
    name = fields.String(required=True)
    identifier = fields.String(required=True)
    namespace = fields.String(required=True)
    measurement = fields.Float(required=True)
    uncertainty = fields.Float(required=True)


class Proteomics(Schema):
    identifier = fields.String(required=True)  # UniProt code
    measurement = fields.Float(required=True)  # mmol/gDw
    uncertainty = fields.Float(required=True)  # mmol/gDw


class UptakeSecretionRates(Schema):
    name = fields.String(required=True)
    identifier = fields.String(required=True)
    namespace = fields.String(required=True)
    measurement = fields.Float(required=True)
    uncertainty = fields.Float(required=True)


class MolarYields(Schema):
    product_name = fields.String(required=True)
    product_identifier = fields.String(required=True)
    product_namespace = fields.String(required=True)
    substrate_name = fields.String(required=True)
    substrate_identifier = fields.String(required=True)
    substrate_namespace = fields.String(required=True)
    measurement = fields.Float(required=True)
    uncertainty = fields.Float(required=True)


class GrowthRate(Schema):
    measurement = fields.Float(required=True)
    uncertainty = fields.Float(required=True)


class ModificationRequest(Schema):
    medium = fields.Nested(MediumCompound, many=True, missing=[])
    genotype = fields.Function(deserialize=gnomic.Genotype.parse, missing=[])
    fluxomics = fields.Nested(Fluxomics, many=True, missing=[])
    metabolomics = fields.Nested(Metabolomics, many=True, missing=[])
    proteomics = fields.Nested(Proteomics, many=True, missing=[])
    uptake_secretion_rates = fields.Nested(UptakeSecretionRates, many=True, missing=[])
    molar_yields = fields.Nested(MolarYields, many=True, missing=[])
    growth_rate = fields.Nested(GrowthRate, missing=None)


class SimulationRequest(Schema):
    model_id = fields.Integer(required=True)
    method = fields.String(missing="fba")
    objective_id = fields.String(missing=None)
    objective_direction = fields.String(missing=None)
    operations = fields.Nested(Operation, many=True, missing=[])
