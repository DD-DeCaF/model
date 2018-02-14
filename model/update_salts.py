import os
import requests
from urllib import request
from multiprocessing import Pool
from functools import partial
import xml.etree.ElementTree as ET


DIR_PATH = os.path.dirname(os.path.realpath(__file__))
FILE_PATH = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/nightly/chebi.obo'
N_PROCESSES = 20


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
        compounds_not_found = missing_compounds(salts, smiles_to_chebi)
        salts = map_salts(smiles_to_chebi)
        print(f'{len(salts)} salts is mapped, '
              f'but one or more compounds are not present for '
              f'{len({k for k, v in salts.items() if [] in v})}')
        # Iron trichloride hack until the proper solution is implemented
        salts[30808] = [[29034], [13291]]
        #
        save_salts(salts)


def save_salts(salts_mapping):
    with open(DIR_PATH + '/data/salts.csv', 'w') as f:
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
