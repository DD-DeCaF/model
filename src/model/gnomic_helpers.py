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
    if not genotype_changes:
        return gnomic.Genotype([])
    genotype = gnomic.Genotype.parse(genotype_changes[0])
    for definition in genotype_changes[1:]:
        genotype = gnomic.Genotype.parse(definition, parent=genotype)
    return genotype


def feature_id(feature):
    """Return the feature identifier (name or accession id) for the given feature"""
    return feature.name or feature.accession.identifier


def feature_additions(genotype_changes):
    """
    Resolve features to add from the given genotype changes

    :param genotype_changes: gnomic.Genotype object with genotype changes
    :return: a generator yielding the names of the features to add
    """
    for change in genotype_changes.changes():
        if isinstance(change, gnomic.Mutation):
            if change.new:
                for feature in change.new.features():
                    yield feature_id(feature)
        if isinstance(change, gnomic.Plasmid):
            for feature in change.features():
                yield feature_id(feature)


def feature_knockouts(genotype_changes):
    """
    Resolve features to knockout from the given genotype changes

    :param genotype_changes: gnomic.Genotype object with genotype changes
    :return: a generator yielding the names of the features to knockout
    """
    for change in genotype_changes.changes():
        if isinstance(change, gnomic.Mutation):
            if change.old:
                for feature in change.old.features():
                    yield feature_id(feature)
