"""
Detects C–H···π interactions between chains by analyzing geometric relationships
between donor C–H groups and aromatic ring planes.

This module:
- Iterates through all candidate C–H donors identified from residue–atom mappings.
- For each donor–ring pair, checks geometric criteria including centroid distance
  and angular orientation (α and β) using ring plane parameters.
- Records all valid C–H···π contacts that satisfy user-defined thresholds
  for distance and angles (C–H–centroid and ring-normal).

Output:
Returns a structured dictionary (`C_H_PI_DICT_STRUCTURE`) storing all detected
C–H···π interactions with residue identifiers, atom names, and geometric values.
"""

import copy

from pi_fucntions import Ring_params
from utils import (
    check_atom_h_pi_condition,
)
from constants import (
    CHAIN1_IDENTIFIER,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    C_H_PI_DIST,
    C_H_PI_THETA1,
    C_H_PI_THETA2,
    C_H_PI_DICT_STRUCTURE,
    C_ATOM_COULUMN,
    C_ATOM_FROM,
    H_ATOM_COULUMN,
    RING_FROM,
    ALPHA_NAME,
    BETA_NAME,
    DISTANCE_NAME
)
from update_constants import tc
C_H_PI_DIST = tc.C_H_PI_DIST
C_H_PI_THETA1 = tc.C_H_PI_THETA1
C_H_PI_THETA2 = tc.C_H_PI_THETA2

def get_c_h_pi_interactions(
    df,
    c_h_combos,
    ring_params: Ring_params,
):
    """Find C–H···π interactions between residues.

    Args:
        df (pd.DataFrame): DataFrame containing residue and atom info.
        c_h_combos (dict): Mapping of chains and residues to C–H atoms.
        ring_params (Ring_params): Ring geometry (centroids, planes).

    Returns:
        dict: All detected C–H···π interactions with geometric details.
    """
    c_h_pi_interactions = copy.deepcopy(C_H_PI_DICT_STRUCTURE)
    # try:
    for chain_c, res_num_c in c_h_combos.items():
        res_num_column_curr = None
        if chain_c == CHAIN1_IDENTIFIER:
            res_num_column_curr = RES_NUM_1
            focus_position = 1
        else:
            res_num_column_curr = RES_NUM_2
            focus_position = 2
        for c_atom_res_rum, all_coords_dicts in res_num_c.items():
            residues_to_check_df = df.loc[df[res_num_column_curr] == c_atom_res_rum] #changed
            for index, row in residues_to_check_df.iterrows():
                for c_atom, coords_dicts in all_coords_dicts.items():
                    c_coords = coords_dicts["c_coords"]
                    for h_atom_dict in coords_dicts["h_atoms"]:
                        h_atom_name = h_atom_dict["name"]
                        h_coords = h_atom_dict["coords"]
                        if focus_position == 1:
                            ring_res_num = row[RES_NUM_2]
                            ring_chain = row[CHAIN_2]
                        else:
                            ring_res_num = row[RES_NUM_1]
                            ring_chain = row[CHAIN_1]
                        if ring_res_num in ring_params.planes_coefficients[ring_chain]:
                            plane_coefficients = ring_params.planes_coefficients[ring_chain][
                                ring_res_num
                            ]
                            centroid = ring_params.centroids[ring_chain][ring_res_num]
                            r_cen, atom_h_centoid_angle, atom_normal_angle = (
                                check_atom_h_pi_condition(
                                    c_coords,
                                    h_coords,
                                    plane_coefficients,
                                    centroid,
                                    C_H_PI_DIST,
                                    C_H_PI_THETA1,
                                    C_H_PI_THETA2,
                                )
                            )
                            if r_cen:
                                if focus_position == 1:
                                    c_h_pi_interactions[RES_NAME_1].append(row[RES_NAME_1])
                                    c_h_pi_interactions[RES_NUM_1].append(row[RES_NUM_1])
                                    c_h_pi_interactions[CHAIN_1].append(row[CHAIN_1])
                                    c_h_pi_interactions[RES_NAME_2].append(row[RES_NAME_2])
                                    c_h_pi_interactions[RES_NUM_2].append(row[RES_NUM_2])
                                    c_h_pi_interactions[CHAIN_2].append(row[CHAIN_2])
                                    c_h_pi_interactions[DISTANCE_NAME].append(round(r_cen,2))
                                    c_h_pi_interactions[C_ATOM_FROM].append(row[RES_NAME_1])
                                    c_h_pi_interactions[RING_FROM].append(
                                        f"{row[RES_NAME_2]}-{row[RES_NUM_2]}"
                                    )
                                else:
                                    c_h_pi_interactions[RES_NAME_1].append(row[RES_NAME_1])
                                    c_h_pi_interactions[RES_NUM_1].append(row[RES_NUM_1])
                                    c_h_pi_interactions[CHAIN_1].append(row[CHAIN_1])
                                    c_h_pi_interactions[RES_NAME_2].append(row[RES_NAME_2])
                                    c_h_pi_interactions[RES_NUM_2].append(row[RES_NUM_2])
                                    c_h_pi_interactions[CHAIN_2].append(row[CHAIN_2])
                                    c_h_pi_interactions[DISTANCE_NAME].append(round(r_cen,2))
                                    c_h_pi_interactions[C_ATOM_FROM].append(row[RES_NAME_2])
                                    c_h_pi_interactions[RING_FROM].append(
                                        f"{row[RES_NAME_1]}-{row[RES_NUM_1]}"
                                    )
                                c_h_pi_interactions[C_ATOM_COULUMN].append(c_atom)
                                c_h_pi_interactions[H_ATOM_COULUMN].append(h_atom_name)
                                c_h_pi_interactions[ALPHA_NAME].append(round(atom_h_centoid_angle,2))
                                c_h_pi_interactions[BETA_NAME].append(round(atom_normal_angle,2))

    return c_h_pi_interactions
    # except Exception as e:
    #     print(f"!!!!!!! could not process because of {e}  !!!!!!\n\n")
    # return
