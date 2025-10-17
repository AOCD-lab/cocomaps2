from identifier_pdb import process_file as process_file_pdb
from identifier_cif import process_file as process_file_cif

def identify_ligands(file_path):
    ligands = {}
    if file_path.endswith(".pdb"):
        ligands = process_file_pdb(file_path)
    elif file_path.endswith(".cif"):
        ligands = process_file_cif(file_path)

    return ligands

