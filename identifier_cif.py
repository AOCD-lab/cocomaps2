import numpy as np

from cif_conversion import Cif_conversion
from read_pdb_regex import read_pdb
from residues import METALS
from constants import ligand_chain, ligand_name, ligand_res_num, rna_chains

LIGAND_EXCLUDES = list(METALS.copy())
LIGAND_EXCLUDES += ["HOH", "BHOH", "AHOH", "3CO", "CL", "2MG",]

threshold = 5


def get_ligands(parsed, ligands):

    ligands_full = []
    for line in parsed:
        if line[3] in ligands and line[0] == "HETATM":
            ligand = (line[3], line[4], line[5])
            if ligand not in ligands_full:
                ligands_full.append(ligand)

    return ligands_full


def process_file(cif_file):
    cif_conv = Cif_conversion(cif_file)
    cif_conv.convert_cif_add_h()
    pdb_file = cif_conv.pdb_file_path
    parsed_lines, raw_lines = read_pdb(
        pdb_file, True, cif_conv.all_res_names, cif_conv.chain_ids
    )
    ligands = []
    for ligand in cif_conv.non_polymors:
        if ligand not in LIGAND_EXCLUDES:
            ligands.append(ligand)
    ligands = get_ligands(parsed_lines, ligands)
    return ligands
# process_file("/workspaces/cocomaps-docker/tests/coco2/5LZS.cif")