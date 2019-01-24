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

import json
import xml.etree.ElementTree as ET
from functools import partial
from multiprocessing import Pool
from urllib import request

import requests


class Chemical:
    """
    chebi: The chebi id, e.g. "CHEBI:12345"
    ions: A list of ions this chemical can be split up to
    ions_missing_smiles: A list of smiles ids that could not be mapped to ions
    metals: A list of metals this chemical can be split up to
    metals_missing_inchi: A list of inchi strings that could not be mapped to metals
    """
    def __init__(self, chebi):
        self.chebi = chebi
        self.ions = set()
        self.ions_missing_smiles = []
        self.metals = set()
        self.metals_missing_inchi = []

    def to_json(self):
        return {
            "ions": list(self.ions),
            "ions_missing_smiles": self.ions_missing_smiles,
            "metals": list(self.metals),
            "metals_missing_inchi": self.metals_missing_inchi,
        }


def main():
    print("Downloading chebi ontology (~125MB)...")
    with request.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo") as file_:
        sections = file_.read().decode("iso-8859-1").split("\n\n")
        print("Done, parsing and creating internal datastructures...")

        # The following will be a complete list of all chemical instances parsed
        # from the chebi ontology
        all_chemicals = []

        # Additionally, create maps for the various identifiers to their
        # corresponding chemical instances
        chebi = {}
        smiles = {}
        inchi = {}
        inchi_formula = {}

        # Parse each "Term" section of the obo file and generate the data
        # structure described above.
        for section in [s for s in sections if s.startswith("[Term]")]:
            parse_obo_term_section(section, all_chemicals, chebi, smiles, inchi, inchi_formula)

    map_ions(all_chemicals, smiles, inchi)
    print(f"  {sum([len(c.ions) for c in all_chemicals])} ions successfully mapped")
    print(f"  {sum([len(c.ions_missing_smiles) for c in all_chemicals])} ions are still unknown")

    resolve_organometallic_metals(all_chemicals, inchi_formula)
    print(f"  {sum([len(c.metals) for c in all_chemicals])} metals mapped")
    print(f"  {sum([len(c.metals_missing_inchi) for c in all_chemicals])} metals are still unknown")

    # TODO (Ali Kaafarani): Deal with nested structures. For example, CHEBI:86368 -> CHEBI:63041 -> CHEBI:29035 etc.
    # Currently, only the topmost level will be mapped to its respective ions/metals.

    # Create a suitable json format and persist it to a file. We're adding indentation and sorting keys
    # for readable diffs at later updates, at a small cost of file size.
    salts = {c.chebi: c.to_json() for c in all_chemicals if c.ions or c.metals}
    with open("data/salts.json", "w") as file_:
        file_.write(json.dumps(salts, indent=2, sort_keys=True))

    print(f"Wrote {len(salts)} salt mappings to 'data/salts.json'")


def parse_obo_term_section(section, all_chemicals, chebi, smiles, inchi, inchi_formula):
    for line in section.split("\n"):
        if line.startswith("id: "):
            # The first line is the identifier; initialize a new chemical instance
            id = line.split()[1]
            chemical = Chemical(id)
            all_chemicals.append(chemical)
            chebi[id] = chemical
        elif line.startswith("property_value"):
            if "smiles" in line:
                # Set the smiles id on the chemical
                id = line.split()[2].strip('"')
                chemical.smiles = id

                # Add the chemical to the smiles map
                if id in smiles:
                    smiles[id].append(chemical)
                else:
                    smiles[id] = [chemical]
            elif "inchi " in line:  # The trailing space separates the inchi string from inchikey
                # Set the inchi values on the chemical
                id = line.split()[2].strip('"')
                chemical.inchi = id
                chemical.inchi_formula = id.split('/', 2)[1]

                # Add the chemical to the inchi formula map
                if chemical.inchi_formula in inchi_formula:
                    inchi_formula[chemical.inchi_formula].append(chemical)
                else:
                    inchi_formula[chemical.inchi_formula] = [chemical]

                # Add the chemical to the inchi map
                if id in inchi:
                    inchi[id].append(chemical)
                else:
                    inchi[id] = [chemical]


def map_ions(all_chemicals, smiles, inchi):
    # Salts are recognized by a period in the smiles id, splitting up the chemicals.
    smiles_salts = [c for c in all_chemicals if hasattr(c, "smiles") and "." in c.smiles]

    # Create a separate map of unique ions, to avoid doing duplicate lookups.
    # Later, we'll map these back to the chemical objects.
    all_ions = {}
    for chemical in smiles_salts:
        for ion in chemical.smiles.split("."):
            all_ions[ion] = None
    print(f"{len(all_ions)} unique ions to map")

    # Map ions using the smiles id
    print("Looking up local smiles maps...")
    for ion in all_ions:
        if ion in smiles:
            all_ions[ion] = smiles[ion]
    print(f"  {len([c for c in all_ions.values() if c is not None])} ions mapped through smiles")

    # Map the remaining ions using inchi keys
    print(f"Mapping chemicals through inchi...")
    print(f"  Looking up {len([i for i, c in all_ions.items() if c is None])} inchi maps in chemspider API (may take a few minutes)...")
    with Pool(processes=20) as pool:
        func = partial(map_inchi, inchi)
        missing_ions = [smiles_ion for smiles_ion, chebi_ids in all_ions.items() if chebi_ids is None]
        for smiles_ion, chemicals in pool.imap_unordered(func, missing_ions, chunksize=100):
            if chemicals is not None:
                all_ions[smiles_ion] = chemicals
    print()

    # Now assign the mapped ions back to the chemical objects
    for chemical in smiles_salts:
        for ion in chemical.smiles.split("."):
            if all_ions[ion] is None:
                chemical.ions_missing_smiles.append(ion)
            else:
                chemical.ions.update([c.chebi for c in all_ions[ion]])


def map_inchi(inchi, smiles_ion):
    """
    Map the given smiles ion to the corresponding inchi string

    Note: The chemspider API lookup should in the future be replaced with local
    conversion using openbabel. See https://pypi.org/project/openbabel/
    """
    print(".", end="", flush=True)
    response = requests.get(f"https://www.chemspider.com/InChI.asmx/SMILESToInChI", params={'smiles': smiles_ion})
    if response.status_code == 500:
        print(f" warning: Could not parse smiles ion: {smiles_ion} ", end="", flush=True)
        return (smiles_ion, None)
    response.raise_for_status()
    inchi_string = ET.fromstring(response.text).text
    try:
        return (smiles_ion, inchi[inchi_string])
    except KeyError:
        return (smiles_ion, None)


def resolve_organometallic_metals(all_chemicals, inchi_formula):
    print("Resolving metals from organometallic compounds...")
    inchi_metals = [c for c in all_chemicals if hasattr(c, "inchi_formula") and "." in c.inchi_formula]
    for chemical in inchi_metals:
        for metal in chemical.inchi_formula.split("."):
            # Check the compound length, stripping the number of elements.
            # We are only looking for metals, so knowing that acid + metal = salt + water,
            # the short compound must be the metal, so we'll look for compounds of 3 chars or
            # less, which might seem a tad arbitrary, but seems to give good results.
            # TODO: The acid should also be added! For example, L-histidine in
            # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32603
            metal_element = metal.lstrip("0123456789")
            if len(metal_element) > 3:
                continue

            if metal_element in inchi_formula:
                chemical.metals.update([c.chebi for c in inchi_formula[metal_element]])
            else:
                chemical.metals_missing_inchi.append(metal_element)


if __name__ == '__main__':
    main()
