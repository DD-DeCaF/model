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

def get_unique_metabolite(model, compound_id, compartment='e', db_name='CHEBI'):
    """Get the only metabolite for given compound / compartment.

    :param model: cobra.Model
    :param compound_id: string, compound identifier, e.g. CHEBI:12965
    :param compartment: string, compartment identifier
    :param db_name: string, the database name, e.g. 'CHEBI'
    """
    # TODO: change id-mapper to use miriam db_names
    key_to_db_name = {'bigg': 'bigg.metabolite', 'chebi': 'CHEBI', 'mnx': 'metanetx.chemical'}
    db_name = key_to_db_name.get(db_name, db_name)
    # TODO: change to only use upper-case chebi everywhere, and always prefix id's with db_name thereby removing the
    # need for the db_name parameter
    if db_name == 'CHEBI':
        compound_id = compound_id.replace('chebi:', 'CHEBI:')
        if re.match('^[0-9]+$', compound_id):
            compound_id = 'CHEBI:' + compound_id

    def query_fun(m):
        xrefs = m.annotation.get(db_name, [])
        xrefs = xrefs if isinstance(xrefs, list) else [xrefs]
        return compound_id in xrefs and m.compartment == compartment

    metabolites = model.metabolites.query(query_fun)
    if len(metabolites) > 1:
        raise IndexError('expected single metabolite, found {}'.format(metabolites))
    if len(metabolites) < 1:
        raise NoIDMapping(compound_id)
    return metabolites[0]


def strip_compartment(x):
    return x[:-2] if re.match('.*_[cepm]$', x) else x
