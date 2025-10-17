"""
Detects weak hydrogen bonds (C–H···O/N type) across interacting chains in biomolecular
complexes. Uses geometric criteria — distance and angle thresholds — to identify
valid donor–acceptor pairs and records them in a structured interaction dictionary.
Builds on atom-level geometry utilities to combine hydrogen bond detection with 
chain-wise parsing and filtering for accurate contact identification.
"""

import copy
from utils import (
    get_acceptor_atom_and_hbond_lines,
    get_angle_between_3_lines,
    get_distance_between_lines,
)

from constants import (
    WEAK_H_BONDS_DICT_STRUCTURE,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_NAME,
    H_ATOM_COULUMN,
    THETA_NAME,
)

from update_constants import tc

CH_ON_DIST = tc.CH_ON_DIST
CH_ON_ANGLE = tc.CH_ON_ANGLE
weak_h_bond_interactions = copy.deepcopy(WEAK_H_BONDS_DICT_STRUCTURE)


def angle_thresholding(line1, line2, h_lines):
    for line3 in h_lines:
        angle = get_angle_between_3_lines(line1, line3, line2)
        if angle >= CH_ON_ANGLE:
            return line3, angle
    return None


def distance_thresholding(line1, line2):
    dist = get_distance_between_lines(line1, line2)
    if dist <= CH_ON_DIST:
        return dist
    return None


def get_weak_h_bond_interactions(chains_parsed_lines, chains):
    weak_h_bond_lines, acceptor_atom_lines = get_acceptor_atom_and_hbond_lines(
        chains_parsed_lines, chains
    )
    for chain, res_weak_h_bond_lines in weak_h_bond_lines.items():
        if chain == chains[0]:
            curr_c_chain_position = 1
            other_chain = chains[1]
        else:
            curr_c_chain_position = 2
            other_chain = chains[0]
        acceptor_atom_lines_other_chain = acceptor_atom_lines[other_chain]
        for res_num, c_collections in res_weak_h_bond_lines.items():
            for c_atom, c_h_lines in c_collections.items():
                for acceptor_atom_line in acceptor_atom_lines_other_chain:
                    satisfying_h_line = angle_thresholding(
                        c_h_lines["C_line"],
                        acceptor_atom_line,
                        c_h_lines["H_lines"],
                    )
                    distance = distance_thresholding(
                        acceptor_atom_line, c_h_lines["C_line"]
                    )
                    if distance and satisfying_h_line:
                        line, angle = satisfying_h_line
                        if curr_c_chain_position == 1:
                            weak_h_bond_interactions[RES_NAME_1].append(
                                c_h_lines["C_line"][3]
                            )
                            weak_h_bond_interactions[RES_NUM_1].append(res_num)
                            weak_h_bond_interactions[CHAIN_1].append(chain)
                            weak_h_bond_interactions[ATOM_1].append(c_atom)
                            weak_h_bond_interactions[RES_NAME_2].append(
                                acceptor_atom_line[3]
                            )
                            weak_h_bond_interactions[RES_NUM_2].append(
                                acceptor_atom_line[5]
                            )
                            weak_h_bond_interactions[CHAIN_2].append(other_chain)
                            weak_h_bond_interactions[ATOM_2].append(
                                acceptor_atom_line[2]
                            )
                        else:
                            weak_h_bond_interactions[RES_NAME_1].append(
                                acceptor_atom_line[3]
                            )
                            weak_h_bond_interactions[RES_NUM_1].append(
                                acceptor_atom_line[5]
                            )
                            weak_h_bond_interactions[CHAIN_1].append(other_chain)
                            weak_h_bond_interactions[ATOM_1].append(
                                acceptor_atom_line[2]
                            )
                            weak_h_bond_interactions[RES_NAME_2].append(
                                c_h_lines["C_line"][3]
                            )
                            weak_h_bond_interactions[RES_NUM_2].append(res_num)
                            weak_h_bond_interactions[CHAIN_2].append(chain)
                            weak_h_bond_interactions[ATOM_2].append(c_atom)
                        weak_h_bond_interactions[DISTANCE_NAME].append(distance)
                        weak_h_bond_interactions[THETA_NAME].append(angle)
                        weak_h_bond_interactions[H_ATOM_COULUMN].append(line[2])
    return weak_h_bond_interactions
