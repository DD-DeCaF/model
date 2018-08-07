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

import gnomic


def full_genotype(genotype_changes):
    """
    Construct gnomic Genotype object from the list of strings with changes

    :param genotype_changes: list of changes, f.e. ['-tyrA::kanMX+', 'kanMX-']
    :return:
    """

    def chain(definitions):
        if not definitions:
            return gnomic.Genotype([])
        genotype = gnomic.Genotype.parse(definitions[0])
        for definition in definitions[1:]:
            genotype = gnomic.Genotype.parse(definition, parent=genotype)
        return genotype

    return chain(genotype_changes)


def insert_feature(feature, dict1, dict2):
    """Helper function for managing two feature storages"""
    if not feature.name:
        return
    if feature.name in dict2:
        dict2.pop(feature.name)
    else:
        dict1[feature.name] = feature


def new_features_identifiers(genotype_changes: gnomic.Genotype):
    """Extract identifiers for features which addition is defined in gnomic string

    :param genotype_changes: gnomic string with genotype changes
    :return:
    """
    for change in genotype_changes.changes():
        if isinstance(change, gnomic.Mutation):
            if change.new:
                for feature in change.new.features():
                    yield feature_id(feature)
        if isinstance(change, gnomic.Plasmid):
            for feature in change.features():
                yield feature_id(feature)
