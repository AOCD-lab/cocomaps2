"""
Detects metal-mediated interactions between residues across chains by
analyzing atomic distances relative to known metal coordination limits.
Uses geometric filtering based on METAL_DISTANCES to identify valid
metal bridges and stores results in a structured dictionary format.
"""

import copy

from utils import (
    get_metal_and_acceptor_lines,
    get_distance_between_lines,
)

from residues import METAL_DISTANCES
from constants import (
    METAL_MEDIATED_DICT_STRUCTURE,
    CHAIN1_IDENTIFIER,
    CHAIN2_IDENTIFIER,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_FROM_RES1_NAME,
    DISTANCE_FROM_RES2_NAME,
    METAL_IDENTITY_NAME,
)
from update_constants import tc


def get_intercations_metal(all_pdb_parsed_lines, chains_mapping, metal_oriented_lines):
    """This fucntion checks the interactions feasibility using METAL_DISTANCES
    dictionary which can be found in the ressidues file.
    Param pdb_file and chains are of same structure as above.
    Return Final_interactions: dict of final thus gotten interactions.
    Structure of final_interactions:

    [metal_atom_line:
                     [donorline1, donorline2],
                     [2,
                     [distance_of metal_line and donorline1, distance_of metal_line and donorline2]
                     ]
                     ]"""
    metal_mediated_interactions = copy.deepcopy(METAL_MEDIATED_DICT_STRUCTURE)
    metal_atoms_lines, acceptor_atom_lines = get_metal_and_acceptor_lines(
        all_pdb_parsed_lines, chains_mapping, metal_oriented_lines
    )
    interactions = []
    for metal_line in metal_atoms_lines:
        # metal_line1 = " ".join(str(e) for e in metal_line)
        metal_identity = f"{metal_line[4]}_{metal_line[5]}_{metal_line[2]}"
        for acceptor_chain1_line in acceptor_atom_lines[CHAIN1_IDENTIFIER]:
            dist1 = get_distance_between_lines(metal_line, acceptor_chain1_line)
            if dist1 <= METAL_DISTANCES[metal_line[11]]:
                for acceptor_chain2_line in acceptor_atom_lines[CHAIN2_IDENTIFIER]:

                    dist2 = get_distance_between_lines(metal_line, acceptor_chain2_line)
                    if dist2 <= METAL_DISTANCES[metal_line[11]]:
                        metal_mediated_interactions[RES_NAME_1].append(
                            acceptor_chain1_line[3]
                        )
                        metal_mediated_interactions[RES_NUM_1].append(
                            acceptor_chain1_line[5]
                        )
                        metal_mediated_interactions[CHAIN_1].append(
                            acceptor_chain1_line[4]
                        )
                        metal_mediated_interactions[ATOM_1].append(
                            acceptor_chain1_line[2]
                        )
                        metal_mediated_interactions[RES_NAME_2].append(
                            acceptor_chain2_line[3]
                        )
                        metal_mediated_interactions[RES_NUM_2].append(
                            acceptor_chain2_line[5]
                        )
                        metal_mediated_interactions[CHAIN_2].append(
                            acceptor_chain2_line[4]
                        )
                        metal_mediated_interactions[ATOM_2].append(
                            acceptor_chain2_line[2]
                        )
                        metal_mediated_interactions[DISTANCE_FROM_RES1_NAME].append(
                            dist1
                        )
                        metal_mediated_interactions[DISTANCE_FROM_RES2_NAME].append(
                            dist2
                        )

                        metal_mediated_interactions[METAL_IDENTITY_NAME].append(
                            metal_identity
                        )

    return metal_mediated_interactions
