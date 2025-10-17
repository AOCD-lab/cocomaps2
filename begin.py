"""
Main execution script for the biomolecular interaction analysis pipeline.

This script serves as the command-line entry point, orchestrating the complete
workflow for detecting inter-chain interactions in protein/RNA complexes.

Workflow overview:
- Parses a user-provided JSON configuration file containing analysis parameters
  and geometric thresholds.
- Optionally converts CIF files to PDB format and adds hydrogens using REDUCE.
- Handles multi-model PDBs by isolating the first model when present.
- Iterates over all chain pair combinations specified in the configuration.
- For each pair, formats PDBs, prepares atom-level data, and runs the main
  `process()` pipeline to compute and save all detected interactions.

Intended usage:
    python main.py config.json

Outputs include detailed CSVs for each interaction type, surface accessibility
tables (via NACCESS), and truncated PDB files for visualization.
"""

import argparse
import time
import os

from Bio.PDB import PDBParser, PDBIO

from update_constants import tc
from cif_conversion import Cif_conversion
from read_pdb_regex import read_pdb
from constants import REDUCE_COMMAND


def parse_arguments():
    """
    Parses command-line arguments and returns them as variables.

    Returns:
        json_file_path (str): Path to the json file with all the threshold values
    """
    parser = argparse.ArgumentParser(description="Process ligand-related information.")

    # Add arguments
    parser.add_argument(
        "json_file_path", type=str, help="Residue number of the ligand (string)"
    )

    args = parser.parse_args()

    return args.json_file_path


if __name__ == "__main__":
    start = time.time()
    json_file_path = parse_arguments()
    tc.read_json_file(json_file_path)

    from workflow import process
    from utils import get_chains_ranges_lines, chain_selection_pdbs_making

    cif_bool = False
    cif_conv = None
    if tc.pdb_file.endswith(".cif"):
        cif_conv = Cif_conversion(tc.pdb_file)
        cif_conv.convert_cif_add_h()
        tc.pdb_file = cif_conv.pdb_file_path
        cif_bool = True

    def check_model_file(pdb_file):
        # new_file = f"{pdb_file[:-4]}_model.pdb"
        fl = open(pdb_file, "r")
        all_file = fl.read()
        if "MODEL        1" in all_file:
            new_fl = open(pdb_file, "w")
            end_ind = all_file.index("ENDMDL")
            model_lines = all_file[:end_ind]
            new_fl.write(model_lines)
            new_fl.close()
            return pdb_file
        else:
            return None

    def format_pdb(pdb_file_path):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("molecule", pdb_file_path)

        # Biopython automatically detects elements from atom names
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file_path)
        return pdb_file_path

    def initial_check(pdb_file_path):
        fl = open(pdb_file_path, "r")
        lines = fl.readlines()
        fl.close()
        for line in lines:
            if line.startswith("ATOM"):
                if line[-1].isalpha():
                    return None
                else:
                    pdb_file_path = format_pdb(pdb_file_path)
                    return pdb_file_path

    new_pdb = initial_check(tc.pdb_file)
    if new_pdb:
        tc.pdb_file = new_pdb
    new_pdb = check_model_file(tc.pdb_file)
    if new_pdb:
        tc.pdb_file = new_pdb

    psrsed, lines = None, None
    if tc.REDUCE_BOOL:
        head, tail = os.path.split(tc.pdb_file)
        new_tail = f"{'.'.join(tail.split('.'))[:-1]}_h.{tail.split('.')[-1]}"
        new_pdb = os.path.join(head, new_tail)
        reduce_command_temp = REDUCE_COMMAND % (tc.pdb_file, new_pdb)
        os.system(reduce_command_temp)
        tc.pdb_file = new_pdb
    if cif_bool:
        parsed, lines = read_pdb(
            tc.pdb_file, cif_bool, cif_conv.all_res_names, cif_conv.chain_ids
        )
    else:
        parsed, lines = read_pdb(tc.pdb_file, cif_bool)

    chain_combinations = []
    for ind1, chain1 in enumerate(tc.chains_set_1):
        for ind2, chain2 in enumerate(tc.chains_set_2):
            if chain1 != chain2:
                try:
                    if (chain1, chain2) not in chain_combinations and (
                        chain2,
                        chain1,
                    ) not in chain_combinations:
                        range_chain_1 = tc.ranges_1[ind1]
                        range_chain_2 = tc.ranges_2[ind2]

                        chain_combinations.append((chain1, chain2))
                        chain_combinations.append((chain2, chain1))
                        (
                            chains_lines,
                            chains_mapping,
                            all_pdb_parsed_lines,
                            all_pdb_formatted_path,
                            compound_pdb_path,
                            chain1_pdb_path,
                            chain2_pdb_path,
                        ) = get_chains_ranges_lines(
                            parsed_lines=parsed,
                            ranges1=range_chain_1,
                            ranges2=range_chain_2,
                            chains=[chain1, chain2],
                            pdb_path=tc.pdb_file,
                        )
                        chain1_pdb_path, chain2_pdb_path, compound_pdb_path = (
                            chain_selection_pdbs_making(chain1, chain2, tc.pdb_file)
                        )

                        process(
                            tc.pdb_file,
                            [chain1, chain2],
                            chains_lines,
                            chains_mapping,
                            all_pdb_parsed_lines,
                            all_pdb_formatted_path,
                            compound_pdb_path,
                            chain1_pdb_path,
                            chain2_pdb_path,
                        )
                except:
                    pass

    end = time.time()
    print(end - start)
