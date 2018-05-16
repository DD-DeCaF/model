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

"""
There could be potential problems with BIGG metabolites that differ in charge
but have the same formula. Keep an eye on:
fe2 and fe3
mn2 and mn4
cu2 and cu
Both from the couple are added to the model, even if only one is the part
of the medium chemical.
"""
import json
import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from urllib import request

import requests


FILE_PATH = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
N_PROCESSES = 20

ID_MAPPER_API = os.environ['ID_MAPPER_API']


def bigg_ids_items(items):
    print("Start querying {} chemicals".format(len(items)))
    return {k: bigg_ids(list(map(str, v))) for k, v in items}


def bigg_ids(object_ids):
    print("Call for {} chemical ids from id mapper".format(len(object_ids)))
    query = json.dumps(
        {'ids': object_ids, 'dbFrom': 'chebi', 'dbTo': 'bigg',
         'type': 'Metabolite'})
    r = requests.post(ID_MAPPER_API, data=query)
    return r.json()['ids']


def create_salts_mapping():
    with request.urlopen(FILE_PATH) as f:
        print('Retrieving the CHEBI file...')
        all_chemicals = [i for i in f.read().decode("iso-8859-1").split('\n\n') if i[:5] == '[Term']
        print('File is ready')
        chebi_to_smiles = identifier_mapping(all_chemicals, 'smiles')
        chebi_to_inchi = identifier_mapping(all_chemicals, 'inchi')
        smiles_to_chebi = revert_mapping(chebi_to_smiles)
        # NH4+ is a part of many media but is not found automatically
        smiles_to_chebi['[NH4+]'] = [28938]
        #
        inchi_to_chebi = revert_mapping(chebi_to_inchi)
        salts = [k for k in smiles_to_chebi if '.' in k]
        compounds_not_found = missing_compounds(salts, smiles_to_chebi)
        with Pool(processes=N_PROCESSES) as pool:
            for result in pool.imap_unordered(
                    partial(smiles_through_inchi, inchi_to_chebi),
                    [compounds_not_found[i::N_PROCESSES] for i in range(N_PROCESSES)]
            ):
                smiles_to_chebi.update(result)
        missing_compounds(salts, smiles_to_chebi)
        salts = map_salts(smiles_to_chebi)
        print(f'{len(salts)} salts is mapped, '
              f'but one or more compounds are not present for '
              f'{len({k for k, v in salts.items() if [] in v})}')

        chebi_to_formula = generate_chebi_to_formula(chebi_to_inchi)
        formula_to_chebi = {}
        for k, v in chebi_to_formula.items():
            formula_to_chebi[v] = formula_to_chebi.get(v, [])
            formula_to_chebi[v].append(k)

        chebi_to_formula_new_salts = {k: v for k, v in chebi_to_formula.items() if '.' in v and k not in salts}

        metal_to_all_chebi = generate_metal_to_chebi(chebi_to_formula_new_salts, formula_to_chebi)
        metal_to_chebi = metal_to_chebi_found_in_bigg(metal_to_all_chebi)

        organometallic_compounds = get_organometallic_compounds(chebi_to_formula_new_salts, metal_to_chebi)

        to_add = defaultdict(list)
        for k, v in salts.items():
            for array in v:
                for element in array:
                    if element in organometallic_compounds:
                        to_add[k].extend(organometallic_compounds[element])
        for k, v in to_add.items():
            salts[k].extend(v)
        salts.update(organometallic_compounds)
        print(f'With {len(organometallic_compounds)} organometallic compounds '
              f'the number of salts is {len(salts)}')

        save_salts(salts)


def get_organometallic_compounds(chebi_to_formula_new_salts,
                                 metal_to_chebi):
    organometallic_compounds = {}
    for k, v in chebi_to_formula_new_salts.items():
        result = []
        for element in v.split('.'):
            element = element.lstrip('0123456789')
            result.append(metal_to_chebi.get(element, []))
        if any(result):
            organometallic_compounds[k] = result
    return organometallic_compounds


def metal_to_chebi_found_in_bigg(metal_to_all_chebi):
    metal_to_bigg = {}
    with Pool(processes=N_PROCESSES) as pool:
        for result in pool.imap_unordered(
                bigg_ids_items,
                [list(metal_to_all_chebi.items())[i::N_PROCESSES] for i
                 in range(N_PROCESSES)]
        ):
            metal_to_bigg.update(result)
    metal_to_chebi = {k: list(map(int, v.keys())) for k, v in
                      metal_to_bigg.items() if v}
    return metal_to_chebi


def generate_metal_to_chebi(chebi_to_formula_new_salts, formula_to_chebi):
    chebi_to_shortest_splits = {}
    all_shortest = []
    for k, v in chebi_to_formula_new_salts.items():
        shortests = [i.lstrip('0123456789') for i in v.split('.') if
                     len(i.lstrip('0123456789')) <= 3]
        all_shortest.extend(shortests)
        chebi_to_shortest_splits[k] = shortests
    all_shortest = set(all_shortest)
    metal_to_all_chebi = {i: formula_to_chebi[i] for i in all_shortest
                          if i in formula_to_chebi}
    return metal_to_all_chebi


def generate_chebi_to_formula(chebi_to_inchi):
    chebi_to_formula = {k: re.match(r'InChI=1S\/([^/]+)\/', v) for k, v in
                        chebi_to_inchi.items()}
    chebi_to_formula = {k: v.group(1) for k, v in chebi_to_formula.items()
                        if v}
    return chebi_to_formula


def save_salts(salts_mapping):
    with open('data/salts.csv', 'w') as f:
        for k, v in salts_mapping.items():
            f.write('{}:{}\n'.format(k, ';'.join(
                [','.join(map(str, i)) for i in v])))


def missing_compounds(salts, smiles_to_chebi):
    print(f'{len(salts)} salts in total')
    compounds_not_found = set()
    compounds_all = set()
    for smiles_string in salts:
        compounds = smiles_string.split('.')
        for comp in compounds:
            compounds_all.add(comp)
            if comp not in smiles_to_chebi:
                compounds_not_found.add(comp)
    print(f'{len(compounds_all)} compounds in total, '
          f'{len(compounds_not_found)} is not found')
    return list(compounds_not_found)


def identifier_mapping(all_chemicals, term):
    mapping = {}
    for row in all_chemicals:
        array = row.split('\n')
        chebi_id = array[1].split(':')[-1]
        entry = [i for i in array if term + ' ' in i]
        if not entry:
            continue
        try:
            key = entry[0].split()[2].strip('"')
        except IndexError:
            print(array)
            raise
        mapping[int(chebi_id)] = key
    return mapping


def revert_mapping(mapping):
    reverse_mapping = {}
    for k, v in mapping.items():
        synonyms = [v, v.replace('+2]', '++]')]
        for syn in synonyms:
            reverse_mapping[syn] = reverse_mapping.get(syn, [])
            if k not in reverse_mapping[syn]:
                reverse_mapping[syn].append(k)
    return reverse_mapping


def smiles_through_inchi(inchi_to_chebi, smiles_strings):
    print(f'Retrieving {len(smiles_strings)} smiles to inchi mappings '
          f'from chemspider')
    smiles_to_chebi = {}
    i = j = 0
    for string in smiles_strings:
        r = requests.get(
            'https://www.chemspider.com/InChI.asmx/SMILESToInChI?smiles={}'.format(string)
        )
        if r.status_code != 500:
            inchi_string = ET.fromstring(r.text).text
            if inchi_string in inchi_to_chebi:
                j += 1
                smiles_to_chebi[string] = inchi_to_chebi[inchi_string]
        i += 1
        if i % 50 == 0:
            print(f'{i} entries queried, {j} found')
    return smiles_to_chebi


def map_salts(smiles_to_chebi):
    salts_mapping = {}
    for smiles_string, chebi_ids in smiles_to_chebi.items():
        if '.' in smiles_string:
            smiles_list = smiles_string.split('.')
            result = []
            for u in smiles_list:
                if u in smiles_to_chebi and smiles_to_chebi[u] not in result:
                    result.append(smiles_to_chebi[u])
                elif u not in smiles_to_chebi:
                    result.append([])
            for ch_id in chebi_ids:
                salts_mapping[ch_id] = result
    return salts_mapping


if __name__ == '__main__':
    create_salts_mapping()
