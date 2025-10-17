"""
Identifies van der Waals (VDW) interactions in molecular complexes, distinguishing
between polar and apolar contacts. Applies atom-type and distance-based criteria
using predefined radii and tolerance thresholds, while excluding atoms already
involved in stronger interactions (e.g., hydrogen bonds, π–π, or halogen bonds).
Returns structured dictionaries of detected polar and apolar VDW contacts.
"""

import copy

from residues import HYDROPHOBIC_NEGLECTS, RADII
from constants import (
    APOLAR_VDW_DICT_STRUCTURE,
    POLAR_VDW_DICT_STRUCTURE,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_1,
    ATOM_2,
    DISTANCE_NAME,
    ATOM_1_ORG_NAME,
    ATOM_2_ORG_NAME
)
from utils import remove_similar_rows_from_dicts
from update_constants import tc

APOLAR_TOLERANCE = tc.APOLAR_TOLERANCE
POLAR_TOLERANCE = tc.POLAR_TOLERANCE


def get_polar_vdw(
    full_df,
    pi_pi_interactions,
    hydrogen_dict,
    halogen_interactions,
    weak_h_bond_interactions,
):
    polar_vdw_dict = copy.deepcopy(POLAR_VDW_DICT_STRUCTURE)
    for index, row in full_df.iterrows():
        try:
            res_1 = row[RES_NAME_1]
            hydrobhobic_atom1 = row[ATOM_1]
            hydrobhobic_atom2 = row[ATOM_2]
            res_num_1 = row[RES_NUM_1]
            chain_1 = row[CHAIN_1]
            res_2 = row[RES_NAME_2]
            res_num_2 = row[RES_NUM_2]
            chain_2 = row[CHAIN_2]
            distance = row[DISTANCE_NAME]
            write_flag = False
            sure_write_flag = False
            if res_1 != "HOH" and res_2 != "HOH":
                if (
                    res_1 in HYDROPHOBIC_NEGLECTS
                    and hydrobhobic_atom1 in HYDROPHOBIC_NEGLECTS[res_1]
                    and hydrobhobic_atom2[0] != "H"
                    and hydrobhobic_atom1[0] != "H"
                    and row[DISTANCE_NAME]
                    <= RADII[row[ATOM_1_ORG_NAME]]
                    + RADII[row[ATOM_2_ORG_NAME]]
                    + POLAR_TOLERANCE
                ):
                    write_flag = True
                elif (
                    # hydrobhobic_atom2[0] == "C"  ##this .......
                    res_2 in HYDROPHOBIC_NEGLECTS
                    and hydrobhobic_atom2 in HYDROPHOBIC_NEGLECTS[res_2]
                    and hydrobhobic_atom2[0] != "H"
                    and hydrobhobic_atom1[0] != "H"
                    and row[DISTANCE_NAME]
                    <= RADII[row[ATOM_1_ORG_NAME]]
                    + RADII[row[ATOM_1_ORG_NAME]]
                    + POLAR_TOLERANCE
                ):
                    write_flag = True
                if (
                    write_flag
                    and res_num_1 not in pi_pi_interactions[RES_NUM_1]
                    and res_num_2 not in pi_pi_interactions[RES_NUM_2]
                ):
                    sure_write_flag = True
                elif res_num_1 in pi_pi_interactions[RES_NUM_1] and write_flag:
                    if (
                        res_num_2
                        != pi_pi_interactions[RES_NUM_2][
                            pi_pi_interactions[RES_NUM_1].index(res_num_1)
                        ]
                    ):
                        sure_write_flag = True

                if (
                    sure_write_flag
                    and hydrobhobic_atom1[0] != "H"
                    and hydrobhobic_atom2[0] != "H"
                ):
                    polar_vdw_dict[RES_NAME_1].append(res_1)
                    polar_vdw_dict[RES_NUM_1].append(res_num_1)
                    polar_vdw_dict[CHAIN_1].append(chain_1)
                    polar_vdw_dict[ATOM_1].append(hydrobhobic_atom1)
                    polar_vdw_dict[RES_NAME_2].append(res_2)
                    polar_vdw_dict[RES_NUM_2].append(res_num_2)
                    polar_vdw_dict[CHAIN_2].append(chain_2)
                    polar_vdw_dict[ATOM_2].append(hydrobhobic_atom2)
                    polar_vdw_dict[DISTANCE_NAME].append(distance)
        except Exception as e:
            print(e)
    polar_vdw_dict = remove_similar_rows_from_dicts(
        polar_vdw_dict, hydrogen_dict, POLAR_VDW_DICT_STRUCTURE
    )
    polar_vdw_dict = remove_similar_rows_from_dicts(
        polar_vdw_dict, halogen_interactions, POLAR_VDW_DICT_STRUCTURE
    )
    polar_vdw_dict = remove_similar_rows_from_dicts(
        polar_vdw_dict, weak_h_bond_interactions, POLAR_VDW_DICT_STRUCTURE
    )

    return polar_vdw_dict


def get_apolar_vdw(full_df, pi_pi_interactions):
    apolar_vdw_dict = copy.deepcopy(APOLAR_VDW_DICT_STRUCTURE)
    for index, row in full_df.iterrows():
        try:
            res_1 = row[RES_NAME_1]
            atom1 = row[ATOM_1]
            atom2 = row[ATOM_2]
            res_num_1 = row[RES_NUM_1]
            chain_1 = row[CHAIN_1]
            res_2 = row[RES_NAME_2]
            res_num_2 = row[RES_NUM_2]
            chain_2 = row[CHAIN_2]
            distance = row[DISTANCE_NAME]
            sure_write_flag = False
            if (
                atom1[0] == "C"  ##this .......
                and atom2[0] == "C"
                and atom1[0] != "H"
                and atom2[0] != "H"
                and res_1 in HYDROPHOBIC_NEGLECTS
                and atom1 not in HYDROPHOBIC_NEGLECTS[res_1]
                and res_2 in HYDROPHOBIC_NEGLECTS
                and atom2 not in HYDROPHOBIC_NEGLECTS[res_2]
                and row[DISTANCE_NAME] <= 2 * RADII["C"] + APOLAR_TOLERANCE
            ):
                if res_num_1 not in pi_pi_interactions[RES_NUM_1]:
                    sure_write_flag = True
                elif res_num_1 in pi_pi_interactions[RES_NUM_1]:
                    if (
                        res_num_2
                        != pi_pi_interactions[RES_NUM_2][
                            pi_pi_interactions[RES_NUM_1].index(res_num_1)
                        ]
                    ):
                        sure_write_flag = True
                if sure_write_flag:
                    apolar_vdw_dict[RES_NAME_1].append(res_1)
                    apolar_vdw_dict[RES_NUM_1].append(res_num_1)
                    apolar_vdw_dict[CHAIN_1].append(chain_1)
                    apolar_vdw_dict[ATOM_1].append(atom1)
                    apolar_vdw_dict[RES_NAME_2].append(res_2)
                    apolar_vdw_dict[RES_NUM_2].append(res_num_2)
                    apolar_vdw_dict[CHAIN_2].append(chain_2)
                    apolar_vdw_dict[ATOM_2].append(atom2)
                    apolar_vdw_dict[DISTANCE_NAME].append(distance)
        except Exception as e:
            print(e)

    return apolar_vdw_dict
