from read_pdb_regex import read_pdb
from residues import METALS
from constants import (
    LIGAND_SEQRES_ACCEPTABLES,
)

LIGAND_EXCLUDES = list(METALS.copy())
LIGAND_EXCLUDES += ["HOH", "BHOH", "AHOH", "3CO", "CL", "2MG",]
threshold = 5


def get_ligands(pdb_file):
    parsed, lines = read_pdb(pdb_file)
    raw_lines = open(pdb_file, "r").readlines()
    ligands = []
    for line in raw_lines:
        if "MODRES" == line.split()[0]:
            LIGAND_EXCLUDES.append(line.split()[2])
        if "HETNAM" == line.split()[0]:
            break
        if "SEQRES" == line.split()[0]:
            parts = line.split()[4:]
            for part in parts:
                if part not in LIGAND_SEQRES_ACCEPTABLES:
                    LIGAND_EXCLUDES.append(part.strip())
    for line in parsed:
        if "HETATM" == line[0]:
            potential_ligand = line[3]
            if potential_ligand not in LIGAND_EXCLUDES:
                ligands.append((potential_ligand, line[4], line[5]))
    ligands = list(set(ligands))
    return set(ligands)

def process_file(pdb_file):
    ligands = get_ligands(pdb_file)
    
    return ligands
# process_file("/workspaces/cocomaps-docker/tests/coco2/3GOT.pdb")