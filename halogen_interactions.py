"""
This module detects and records halogen bonding interactions between residues
in molecular structures. It identifies potential halogen donors and acceptors
from parsed chain data, calculates geometric parameters such as distances and
angles, and validates interactions based on predefined thresholds for distance
and angle criteria. The results are returned in a structured dictionary format.
"""

import copy

from utils import (
    get_residue_chunks,
    get_closest_atom_in_residue,
    get_angle_between_3_points,
    get_distance_between_lines,
)
from update_constants import tc
from constants import (
    HALOGEN_DICT_STRUCTURE,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_NAME,
    THETA_1_NAME,
    THETA_2_NAME,
    CHAIN1_IDENTIFIER,
)
from residues import HALOGEN_ELECTRONEGATIVE, HALOGENS, RADII

HALOGEN_THETA1 = tc.HALOGEN_THETA1
HALOGEN_THETA2 = tc.HALOGEN_THETA2


def get_focus_lines(chains_parsed_lines, chains):
    halogen_lines = {chains[0]: [], chains[1]: []}
    acceptor_lines = {chains[0]: [], chains[1]: []}
    residue_chunks = get_residue_chunks(chains_parsed_lines, chains)
    for line in chains_parsed_lines:
        if line[4] in chains:
            if line[11] in HALOGENS:
                halogen_lines[line[4]].append(line)
            if (
                line[4] in chains
                and line[3] in HALOGEN_ELECTRONEGATIVE
                and line[2] in HALOGEN_ELECTRONEGATIVE[line[3]]
            ):
                acceptor_lines[line[4]].append(line)
    return halogen_lines, acceptor_lines, residue_chunks


def check_interactions_feasibility(
    halogen_lines,
    acceptor_lines,
    residue_chunks,
    halogen_dict,
):
    for curr_halogen_chain, curr_chain_halogen_lines in halogen_lines.items():
        if curr_halogen_chain == CHAIN1_IDENTIFIER:
            halogen_res_1_bool = True
        else:
            halogen_res_1_bool = False
        for halogen_line in curr_chain_halogen_lines:
            halogen_closest_atom = None
            halogen_closest_atom = get_closest_atom_in_residue(
                halogen_line, residue_chunks[curr_halogen_chain][halogen_line[5]]
            )
            for (
                curr_acceptor_chain,
                curr_chain_acceptor_lines,
            ) in acceptor_lines.items():
                if curr_acceptor_chain != curr_halogen_chain and halogen_closest_atom:
                    for acceptor_line in curr_chain_acceptor_lines:
                        acceptor_closest_atom = None
                        acceptor_closest_atom = get_closest_atom_in_residue(
                            acceptor_line,
                            residue_chunks[curr_acceptor_chain][acceptor_line[5]],
                        )
                        acceptor_coords = acceptor_line[6:9]
                        halogen_coords = halogen_line[6:9]
                        halogen_closest_coords = halogen_closest_atom[6:9]
                        acceptor_closest_coords = acceptor_closest_atom[6:9]
                        halogen_theta1 = get_angle_between_3_points(
                            acceptor_coords, halogen_coords, halogen_closest_coords
                        )
                        halogen_theta2 = get_angle_between_3_points(
                            halogen_coords, acceptor_coords, acceptor_closest_coords
                        )
                        distance = get_distance_between_lines(
                            halogen_line, acceptor_line
                        )
                        atom1 = "".join(filter(str.isalpha, halogen_line[11]))
                        atom2 = "".join(filter(str.isalpha, acceptor_line[11]))
                        if (
                            halogen_theta1 >= HALOGEN_THETA1
                            and halogen_theta2 >= HALOGEN_THETA2
                            and distance <= RADII[atom1] + RADII[atom2]
                        ):
                            if not halogen_res_1_bool:
                                halogen_dict[RES_NAME_1].append(acceptor_line[3])
                                halogen_dict[RES_NUM_1].append(acceptor_line[5])
                                halogen_dict[CHAIN_1].append(acceptor_line[4])
                                halogen_dict[ATOM_1].append(acceptor_line[2])
                                halogen_dict[RES_NAME_2].append(halogen_line[3])
                                halogen_dict[RES_NUM_2].append(halogen_line[5])
                                halogen_dict[CHAIN_2].append(halogen_line[4])
                                halogen_dict[ATOM_2].append(halogen_line[2])
                                halogen_dict[DISTANCE_NAME].append(distance)
                                halogen_dict[THETA_1_NAME].append(halogen_theta1)
                                halogen_dict[THETA_2_NAME].append(halogen_theta2)
                            else:
                                halogen_dict[RES_NAME_1].append(halogen_line[3])
                                halogen_dict[RES_NUM_1].append(halogen_line[5])
                                halogen_dict[CHAIN_1].append(halogen_line[4])
                                halogen_dict[ATOM_1].append(halogen_line[2])
                                halogen_dict[RES_NAME_2].append(acceptor_line[3])
                                halogen_dict[RES_NUM_2].append(acceptor_line[5])
                                halogen_dict[CHAIN_2].append(acceptor_line[4])
                                halogen_dict[ATOM_2].append(acceptor_line[2])
                                halogen_dict[DISTANCE_NAME].append(distance)
                                halogen_dict[THETA_1_NAME].append(halogen_theta1)
                                halogen_dict[THETA_2_NAME].append(halogen_theta2)
    return halogen_dict


def get_halogen_interactions(chains_parsed_lines, chains):
    halogen_dict = copy.deepcopy(HALOGEN_DICT_STRUCTURE)
    (
        halogen_lines,
        acceptor_lines,
        residue_chunks,
    ) = get_focus_lines(chains_parsed_lines, chains)
    halogen_dict = check_interactions_feasibility(
        halogen_lines,
        acceptor_lines,
        residue_chunks,
        halogen_dict,
    )
    return halogen_dict
