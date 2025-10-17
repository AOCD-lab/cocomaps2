"""
Identifies N/S/O–H···π interactions between donor groups (N, S, or O with attached H)
and aromatic ring systems in molecular structures. Calculates geometric parameters
such as distances and angles using ring plane data to validate interactions, and
returns them in a structured dictionary format.
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
    NSOH_PI_DIST,
    NSOH_PI_THETA1,
    NSOH_PI_THETA2,
    NSO_H_PI_DICT_STRUCTURE,
    C_ATOM_COULUMN,
    C_ATOM_FROM,
    H_ATOM_COULUMN,
    RING_FROM,
    ALPHA_NAME,
    BETA_NAME,
    DISTANCE_NAME
)
from update_constants import tc
NSOH_PI_DIST = tc.NSOH_PI_DIST
NSOH_PI_THETA1 = tc.NSOH_PI_THETA1
NSOH_PI_THETA2 = tc.NSOH_PI_THETA2


def get_nso_h_pi_interactions(
    df,
    nso_h_combos,
    ring_params: Ring_params,
):
    nso_h_pi_interactions = copy.deepcopy(NSO_H_PI_DICT_STRUCTURE)
    # try:
    for chain_nso, res_num_nso in nso_h_combos.items():
        res_num_column_curr = None
        if chain_nso == CHAIN1_IDENTIFIER:
            res_num_column_curr = RES_NUM_1
            focus_position = 1
        else:
            res_num_column_curr = RES_NUM_2
            focus_position = 2
        for nso_atom_res_num, all_coords_dicts in res_num_nso.items():
            residues_to_check_df = df.loc[df[res_num_column_curr] == nso_atom_res_num]
            for index, row in residues_to_check_df.iterrows():
                for nso_atom, coords_dicts in all_coords_dicts.items():
                    nso_coords = coords_dicts["n_coords"]
                    for h_atom_dict in coords_dicts["h_atoms"]:
                        h_atom_name = h_atom_dict["name"]
                        h_coords = h_atom_dict["coords"]
                        if focus_position == 1:
                            ring_res_num = row[RES_NUM_2]
                            ring_chain = row[CHAIN_2]
                        else:
                            ring_res_num = row[RES_NUM_1]
                            ring_chain = row[CHAIN_1]
                        if ring_res_num in plane_coefficients[ring_chain]:
                            plane_coefficients = ring_params.planes_coefficients[ring_chain][
                                ring_res_num
                            ]
                            centroid = ring_params.centroids[ring_chain][ring_res_num]
                            r_cen, atom_h_centoid_angle, atom_normal_angle = (
                                check_atom_h_pi_condition(
                                    nso_coords,
                                    h_coords,
                                    plane_coefficients,
                                    centroid,
                                    NSOH_PI_DIST,
                                    NSOH_PI_THETA1,
                                    NSOH_PI_THETA2,
                                )
                            )
                            if r_cen:
                                if focus_position == 1:
                                    nso_h_pi_interactions[RES_NAME_1].append(row[RES_NAME_1])
                                    nso_h_pi_interactions[RES_NUM_1].append(row[RES_NUM_1])
                                    nso_h_pi_interactions[CHAIN_1].append(row[CHAIN_1])
                                    nso_h_pi_interactions[RES_NAME_2].append(row[RES_NAME_2])
                                    nso_h_pi_interactions[RES_NUM_2].append(row[RES_NUM_2])
                                    nso_h_pi_interactions[CHAIN_2].append(row[CHAIN_2])
                                    nso_h_pi_interactions[DISTANCE_NAME].append(round(r_cen,2))
                                    nso_h_pi_interactions[C_ATOM_FROM].append(row[RES_NAME_1])
                                    nso_h_pi_interactions[RING_FROM].append(
                                        f"{row[RES_NAME_2]}-{row[RES_NUM_2]}"
                                    )
                                else:
                                    nso_h_pi_interactions[RES_NAME_1].append(row[RES_NAME_1])
                                    nso_h_pi_interactions[RES_NUM_1].append(row[RES_NUM_1])
                                    nso_h_pi_interactions[CHAIN_1].append(row[CHAIN_1])
                                    nso_h_pi_interactions[RES_NAME_2].append(row[RES_NAME_2])
                                    nso_h_pi_interactions[RES_NUM_2].append(row[RES_NUM_2])
                                    nso_h_pi_interactions[CHAIN_2].append(row[CHAIN_2])
                                    nso_h_pi_interactions[DISTANCE_NAME].append(round(r_cen,2))
                                    nso_h_pi_interactions[C_ATOM_FROM].append(row[RES_NAME_2])
                                    nso_h_pi_interactions[RING_FROM].append(
                                        f"{row[RES_NAME_1]}-{row[RES_NUM_1]}"
                                    )
                                nso_h_pi_interactions[C_ATOM_COULUMN].append(nso_atom)
                                nso_h_pi_interactions[H_ATOM_COULUMN].append(h_atom_name)
                                nso_h_pi_interactions[ALPHA_NAME].append(round(atom_h_centoid_angle,2))
                                nso_h_pi_interactions[BETA_NAME].append(round(atom_normal_angle,2))

    return nso_h_pi_interactions
    # except Exception as e:
    #     print(f"!!!!!!! could not process because of {e}  !!!!!!\n\n")
