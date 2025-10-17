"""
This module handles conversion and preprocessing of CIF (Crystallographic Information File) structures.
It extracts residue and chain information, converts CIF files into PDB format, and can optionally
add hydrogens using the Reduce tool. Useful for preparing structural data for further molecular analysis.
"""

import os

from Bio.PDB import MMCIFParser, PDBParser
from iotbx.command_line import cif_as_pdb

from read_pdb_regex import read_pdb

from constants import REDUCE_COMMAND


class Cif_conversion:
    """Handles conversion of CIF files to PDB, extraction of residue info, and hydrogen addition."""
    def __init__(self, cif_file_path) -> None:
        self.cif_file_path = cif_file_path
        self.non_polymors = []

    def get_all_res_names_from_file(self):
        """Extract all residue names and non-polymer entities from the CIF file."""
        fl = open(self.cif_file_path)
        lines = fl.readlines()
        all_res_names = []
        flag_names_start = False
        flag_names_start1 = False
        flag_names_start2 = False
        for index, line in enumerate(
            lines
        ):  # added new place to get res names from to get 1PE for 3dil
            if "_entity_poly_seq.hetero" in line:
                flag_names_start = True
                continue
            if flag_names_start and "#" not in line:
                all_res_names.append(line.split()[2].strip())
            elif flag_names_start and "#" in line:
                flag_names_start = False

            if "_chem_comp.id" in line:
                flag_names_start2 = True
                continue
            if flag_names_start2 and "#" not in line and "non-polymer" in line:
                entity = line.split()[0].strip()
                all_res_names.append(entity)
                self.non_polymors.append(entity)
            elif flag_names_start2 and "#" in line:
                flag_names_start2 = False

            if "_pdbx_entity_nonpoly.comp_id" in line:
                flag_names_start1 = True
                continue
            if flag_names_start1 and "#" not in line:
                all_res_names.append(line.split()[-1].strip())
            elif flag_names_start1 and "#" in line:
                flag_names_start1 = False
            if "ATOM " in line:
                break
        self.non_polymors = list(set(self.non_polymors))
        self.all_res_names = list(set(all_res_names))

    def get_all_rna_chains(self):
        """Parse CIF file to extract all chain identifiers."""
        pdb_name = os.path.basename(self.cif_file_path)
        pdb = MMCIFParser().get_structure(pdb_name, self.cif_file_path)
        chain_ids = []
        for chain in pdb.get_chains():
            chain_ids.append(chain.get_id())
        self.chain_ids = chain_ids

    def convert_cif_to_pdb(self):
        self.pdb_file_path = f"{self.cif_file_path[:-4]}.pdb"
        cif_as_pdb.run([self.cif_file_path], self.pdb_file_path)

    def run_reduce(self):
        head, tail = os.path.split(self.pdb_file_path)
        new_tail = f"{'.'.join(tail.split('.'))[:-1]}_h.{tail.split('.')[-1]}"
        new_pdb = os.path.join(head, new_tail)
        os.system(REDUCE_COMMAND % (self.pdb_file_path, new_pdb))
        pdb_file = new_pdb
        self.hydrogen_added_pdb = pdb_file

    def convert_cif_add_h(self):
        """Full pipeline: extract info, convert to PDB, optionally add hydrogens."""
        self.get_all_res_names_from_file()
        self.get_all_rna_chains()
        self.convert_cif_to_pdb()
        # self.run_reduce()
        # return self.hydrogen_added_pdb
