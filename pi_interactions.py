"""
Defines the Ring_params class, which computes and stores geometric parameters for
aromatic ring systems detected in molecular structures. The class calculates each
ring’s centroid, plane coefficients, and polygonal boundaries based on atom coordinates,
allowing accurate geometric reasoning about π-related interactions.

It also provides methods to project nearby atoms or points onto ring planes and
evaluate whether those points lie within or near the π-system boundary, based on
distance thresholds specific to interaction types such as cation–π, anion–π, or
lone pair–π contacts. This enables the geometric validation of π-centered
noncovalent interactions in complex biomolecular assemblies.
"""

import numpy as np
import copy
import math

from pi_fucntions import Ring_params
from utils import get_dist_between_points
from constants import (
    LONE_PAIR_PI_DICT_STRUCTURE,
    ANION_PI_DICT_STRUCTURE,
    CATION_PI_DICT_STRUCTURE,
    AMINO_PI_DICT_STRUCTURE,
    PI_PI_DICT_STRUCTURE,
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
    DISTANCE_NAME,
    THETA_NAME,
    GAMA_NAME,
    RING_FROM,
    LONE_PAIR_ATOM,
    LONE_PAIR_FROM,
    ANION_ATOM,
    ANION_FROM,
    CATION_ATOM,
    CATION_FROM,
    POLAR_ATOM,
    POLAR_FROM,
    LONE_PAIR_PI_PRINT_NAME,
    ANION_PI_PRINT_NAME,
    CATION_PI_PRINT_NAME,
    AMINO_PI_PRINT_NAME,
)
from update_constants import tc
from residues import lone_pair_pi, anion_pi, cation_pi, amino_pi, residues

PI_PI_DIST = tc.PI_PI_DIST
PI_PI_THETA = tc.PI_PI_THETA
PI_PI_GAMMA = tc.PI_PI_GAMMA

lone_pair_pi_interactions = copy.deepcopy(LONE_PAIR_PI_DICT_STRUCTURE)
anion_pi_interactions = copy.deepcopy(ANION_PI_DICT_STRUCTURE)
cation_pi_interactions = copy.deepcopy(CATION_PI_DICT_STRUCTURE)
amino_pi_interactions = copy.deepcopy(AMINO_PI_DICT_STRUCTURE)


def get_pi_pi_interactions(
    ring_params: Ring_params,
    df,
):
    pi_pi_interactions = copy.deepcopy(PI_PI_DICT_STRUCTURE)
    chain1_res_num_dict = ring_params.centroids[CHAIN1_IDENTIFIER]
    chain2_res_num_dict = ring_params.centroids[CHAIN2_IDENTIFIER]
    for res_num1, ring_centroid1 in chain1_res_num_dict.items():
        for res_num2, ring_centroid2 in chain2_res_num_dict.items():
            r_cen = get_dist_between_points(ring_centroid1, ring_centroid2)
            if r_cen <= PI_PI_DIST:
                plane1 = ring_params.planes_coefficients[CHAIN1_IDENTIFIER][res_num1]
                plane2 = ring_params.planes_coefficients[CHAIN2_IDENTIFIER][res_num2]
                normal1 = plane1[0], plane1[1], plane1[2]
                normal2 = plane2[0], plane2[1], plane2[2]
                centroid_vector = np.array(ring_centroid1) - np.array(ring_centroid2)
                dot_prod = np.linalg.norm(np.dot(normal1, centroid_vector))
                cos_theta = dot_prod / (
                    r_cen
                    * math.sqrt(
                        plane1[0] * plane1[0]
                        + plane1[1] * plane1[1]
                        + plane1[2] * plane1[2]
                    )
                )
                theta = np.degrees(np.arccos(cos_theta))
                dot_prod_ = np.linalg.norm(np.dot(normal2, normal1))
                cos_theta_ = dot_prod_ / (
                    r_cen
                    * math.sqrt(
                        plane2[0] * plane2[0]
                        + plane2[1] * plane2[1]
                        + plane2[2] * plane2[2]
                    )
                )
                gama = np.degrees(np.arccos(cos_theta_))
                if theta <= PI_PI_THETA and gama <= PI_PI_GAMMA:
                    if df.loc[
                            (df[RES_NUM_1] == res_num1)
                            & (df[CHAIN_1] == CHAIN1_IDENTIFIER),
                            RES_NAME_1,
                        ].any():
                        res_name1 = df.loc[
                            (df[RES_NUM_1] == res_num1)
                            & (df[CHAIN_1] == CHAIN1_IDENTIFIER),
                            RES_NAME_1,
                        ].values[0]
                        pi_pi_interactions[RES_NAME_1].append(
                            res_name1
                        )
                        pi_pi_interactions[RES_NUM_1].append(res_num1)
                        pi_pi_interactions[CHAIN_1].append(CHAIN1_IDENTIFIER)
                        pi_pi_interactions[RES_NAME_2].append(
                            df.loc[
                                (df[RES_NUM_2] == res_num2)
                                & (df[CHAIN_2] == CHAIN2_IDENTIFIER),
                                RES_NAME_2,
                            ].values[0]
                        )
                        pi_pi_interactions[RES_NUM_2].append(res_num2)
                        pi_pi_interactions[CHAIN_2].append(CHAIN2_IDENTIFIER)
                        pi_pi_interactions[DISTANCE_NAME].append(r_cen)
                        pi_pi_interactions[THETA_NAME].append(theta)
                        pi_pi_interactions[GAMA_NAME].append(gama)
    return pi_pi_interactions


def check_lone_pair_pi(
    atom_coords, atom_row, ring_row, ring_params: Ring_params, chain1_atom_bool: bool
):
    if chain1_atom_bool:
        ring_chain = ring_row[CHAIN_2]
        ring_res_num = ring_row[RES_NUM_2]
    else:
        ring_chain = ring_row[CHAIN_1]
        ring_res_num = ring_row[RES_NUM_1]

    if (
        ring_res_num not in ring_params.planes_coefficients[ring_chain]
    ):  # check with bhaiya
        return

    plane_coefficients = ring_params.planes_coefficients[ring_chain][ring_res_num]
    polygon_coords = ring_params.polygons_coords[ring_chain][ring_res_num]
    r_cen = ring_params.check_point_inside(
        atom_coords, plane_coefficients, polygon_coords, LONE_PAIR_PI_PRINT_NAME, ring_params.centroids[ring_chain][ring_res_num]
    )
    if r_cen and chain1_atom_bool:
        lone_pair_pi_interactions[RES_NAME_1].append(atom_row[RES_NAME_1])
        lone_pair_pi_interactions[RES_NUM_1].append(atom_row[RES_NUM_1])
        lone_pair_pi_interactions[CHAIN_1].append(atom_row[CHAIN_1])
        lone_pair_pi_interactions[RES_NAME_2].append(ring_row[RES_NAME_2])
        lone_pair_pi_interactions[RES_NUM_2].append(ring_row[RES_NUM_2])
        lone_pair_pi_interactions[CHAIN_2].append(ring_row[CHAIN_2])
        lone_pair_pi_interactions[DISTANCE_NAME].append(r_cen)
        lone_pair_pi_interactions[LONE_PAIR_ATOM].append(atom_row[ATOM_1])
        lone_pair_pi_interactions[LONE_PAIR_FROM].append(
            f"{atom_row[RES_NAME_1]}-{atom_row[RES_NUM_1]}"
        )
        lone_pair_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_2]}-{ring_row[RES_NUM_2]}"
        )
    elif r_cen and chain1_atom_bool == False:
        lone_pair_pi_interactions[RES_NAME_1].append(ring_row[RES_NAME_1])
        lone_pair_pi_interactions[RES_NUM_1].append(ring_row[RES_NUM_1])
        lone_pair_pi_interactions[CHAIN_1].append(ring_row[CHAIN_1])
        lone_pair_pi_interactions[RES_NAME_2].append(atom_row[RES_NAME_2])
        lone_pair_pi_interactions[RES_NUM_2].append(atom_row[RES_NUM_2])
        lone_pair_pi_interactions[CHAIN_2].append(atom_row[CHAIN_2])
        lone_pair_pi_interactions[DISTANCE_NAME].append(r_cen)
        lone_pair_pi_interactions[LONE_PAIR_ATOM].append(atom_row[ATOM_2])
        lone_pair_pi_interactions[LONE_PAIR_FROM].append(
            f"{atom_row[RES_NAME_2]}-{atom_row[RES_NUM_2]}"
        )
        lone_pair_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_1]}-{ring_row[RES_NUM_1]}"
        )

    pass


def check_anion_pi(
    atom_coords, atom_row, ring_row, ring_params: Ring_params, chain1_atom_bool: bool
):
    if chain1_atom_bool:
        ring_chain = ring_row[CHAIN_2]
        ring_res_num = ring_row[RES_NUM_2]
    else:
        ring_chain = ring_row[CHAIN_1]
        ring_res_num = ring_row[RES_NUM_1]

    if ring_res_num not in ring_params.planes_coefficients[ring_chain]:
        return

    plane_coefficients = ring_params.planes_coefficients[ring_chain][ring_res_num]
    polygon_coords = ring_params.polygons_coords[ring_chain][ring_res_num]
    r_cen = ring_params.check_point_inside(
        atom_coords, plane_coefficients, polygon_coords, ANION_PI_PRINT_NAME, ring_params.centroids[ring_chain][ring_res_num]
    )
    if r_cen and chain1_atom_bool:
        anion_pi_interactions[RES_NAME_1].append(atom_row[RES_NAME_1])
        anion_pi_interactions[RES_NUM_1].append(atom_row[RES_NUM_1])
        anion_pi_interactions[CHAIN_1].append(atom_row[CHAIN_1])
        anion_pi_interactions[RES_NAME_2].append(ring_row[RES_NAME_2])
        anion_pi_interactions[RES_NUM_2].append(ring_row[RES_NUM_2])
        anion_pi_interactions[CHAIN_2].append(ring_row[CHAIN_2])
        anion_pi_interactions[DISTANCE_NAME].append(r_cen)
        anion_pi_interactions[ANION_ATOM].append(atom_row[ATOM_1])
        anion_pi_interactions[ANION_FROM].append(
            f"{atom_row[RES_NAME_1]}-{atom_row[RES_NUM_1]}"
        )
        anion_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_2]}-{ring_row[RES_NUM_2]}"
        )
    elif r_cen and chain1_atom_bool == False:
        anion_pi_interactions[RES_NAME_1].append(ring_row[RES_NAME_1])
        anion_pi_interactions[RES_NUM_1].append(ring_row[RES_NUM_1])
        anion_pi_interactions[CHAIN_1].append(ring_row[CHAIN_1])
        anion_pi_interactions[RES_NAME_2].append(atom_row[RES_NAME_2])
        anion_pi_interactions[RES_NUM_2].append(atom_row[RES_NUM_2])
        anion_pi_interactions[CHAIN_2].append(atom_row[CHAIN_2])
        anion_pi_interactions[DISTANCE_NAME].append(r_cen)
        anion_pi_interactions[ANION_ATOM].append(atom_row[ATOM_2])
        anion_pi_interactions[ANION_FROM].append(
            f"{atom_row[RES_NAME_2]}-{atom_row[RES_NUM_2]}"
        )
        anion_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_1]}-{ring_row[RES_NUM_1]}"
        )


def check_cation_pi(
    atom_coords, atom_row, ring_row, ring_params: Ring_params, chain1_atom_bool: bool
):
    if chain1_atom_bool:
        ring_chain = ring_row[CHAIN_2]
        ring_res_num = ring_row[RES_NUM_2]
    else:
        ring_chain = ring_row[CHAIN_1]
        ring_res_num = ring_row[RES_NUM_1]
    if ring_res_num not in ring_params.planes_coefficients[ring_chain]:
        return
    plane_coefficients = ring_params.planes_coefficients[ring_chain][ring_res_num]
    polygon_coords = ring_params.polygons_coords[ring_chain][ring_res_num]
    r_cen = ring_params.check_point_inside(
        atom_coords, plane_coefficients, polygon_coords, CATION_PI_PRINT_NAME, ring_params.centroids[ring_chain][ring_res_num]
    )
    if r_cen and chain1_atom_bool:
        cation_pi_interactions[RES_NAME_1].append(atom_row[RES_NAME_1])
        cation_pi_interactions[RES_NUM_1].append(atom_row[RES_NUM_1])
        cation_pi_interactions[CHAIN_1].append(atom_row[CHAIN_1])
        cation_pi_interactions[RES_NAME_2].append(ring_row[RES_NAME_2])
        cation_pi_interactions[RES_NUM_2].append(ring_row[RES_NUM_2])
        cation_pi_interactions[CHAIN_2].append(ring_row[CHAIN_2])
        cation_pi_interactions[DISTANCE_NAME].append(r_cen)
        cation_pi_interactions[CATION_ATOM].append(atom_row[ATOM_1])
        cation_pi_interactions[CATION_FROM].append(
            f"{atom_row[RES_NAME_1]}-{atom_row[RES_NUM_1]}"
        )
        cation_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_2]}-{ring_row[RES_NUM_2]}"
        )
    elif r_cen and chain1_atom_bool == False:
        cation_pi_interactions[RES_NAME_1].append(ring_row[RES_NAME_1])
        cation_pi_interactions[RES_NUM_1].append(ring_row[RES_NUM_1])
        cation_pi_interactions[CHAIN_1].append(ring_row[CHAIN_1])
        cation_pi_interactions[RES_NAME_2].append(atom_row[RES_NAME_2])
        cation_pi_interactions[RES_NUM_2].append(atom_row[RES_NUM_2])
        cation_pi_interactions[CHAIN_2].append(atom_row[CHAIN_2])
        cation_pi_interactions[DISTANCE_NAME].append(r_cen)
        cation_pi_interactions[CATION_ATOM].append(atom_row[ATOM_2])
        cation_pi_interactions[CATION_FROM].append(
            f"{atom_row[RES_NAME_2]}-{atom_row[RES_NUM_2]}"
        )
        cation_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_1]}-{ring_row[RES_NUM_1]}"
        )


def check_amino_pi(
    atom_coords, atom_row, ring_row, ring_params: Ring_params, chain1_atom_bool: bool
):
    if chain1_atom_bool:
        ring_chain = ring_row[CHAIN_2]
        ring_res_num = ring_row[RES_NUM_2]
    else:
        ring_chain = ring_row[CHAIN_1]
        ring_res_num = ring_row[RES_NUM_1]
    if (
        ring_res_num not in ring_params.planes_coefficients[ring_chain]
    ):  # check with bhaiya
        return
    plane_coefficients = ring_params.planes_coefficients[ring_chain][ring_res_num]
    polygon_coords = ring_params.polygons_coords[ring_chain][ring_res_num]
    r_cen = ring_params.check_point_inside(
        atom_coords, plane_coefficients, polygon_coords, AMINO_PI_PRINT_NAME, ring_params.centroids[ring_chain][ring_res_num]
    )
    if r_cen and chain1_atom_bool:
        amino_pi_interactions[RES_NAME_1].append(atom_row[RES_NAME_1])
        amino_pi_interactions[RES_NUM_1].append(atom_row[RES_NUM_1])
        amino_pi_interactions[CHAIN_1].append(atom_row[CHAIN_1])
        amino_pi_interactions[RES_NAME_2].append(ring_row[RES_NAME_2])
        amino_pi_interactions[RES_NUM_2].append(ring_row[RES_NUM_2])
        amino_pi_interactions[CHAIN_2].append(ring_row[CHAIN_2])
        amino_pi_interactions[DISTANCE_NAME].append(r_cen)
        amino_pi_interactions[POLAR_ATOM].append(atom_row[ATOM_1])
        amino_pi_interactions[POLAR_FROM].append(
            f"{atom_row[RES_NAME_1]}-{atom_row[RES_NUM_1]}"
        )
        amino_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_2]}-{ring_row[RES_NUM_2]}"
        )
    elif r_cen and chain1_atom_bool == False:
        amino_pi_interactions[RES_NAME_1].append(ring_row[RES_NAME_1])
        amino_pi_interactions[RES_NUM_1].append(ring_row[RES_NUM_1])
        amino_pi_interactions[CHAIN_1].append(ring_row[CHAIN_1])
        amino_pi_interactions[RES_NAME_2].append(atom_row[RES_NAME_2])
        amino_pi_interactions[RES_NUM_2].append(atom_row[RES_NUM_2])
        amino_pi_interactions[CHAIN_2].append(atom_row[CHAIN_2])
        amino_pi_interactions[DISTANCE_NAME].append(r_cen)
        amino_pi_interactions[POLAR_ATOM].append(atom_row[ATOM_2])
        amino_pi_interactions[POLAR_FROM].append(
            f"{atom_row[RES_NAME_2]}-{atom_row[RES_NUM_2]}"
        )
        amino_pi_interactions[RING_FROM].append(
            f"{ring_row[RES_NAME_1]}-{ring_row[RES_NUM_1]}"
        )


def get_all_pi_interactions(full_df, ring_params: Ring_params, df):

    for index, row in full_df.iterrows():
        res_name1 = row[RES_NAME_1]
        res_num1 = row[RES_NUM_1]
        chain1 = row[CHAIN_1]
        atom1 = row[ATOM_1]
        res_name2 = row[RES_NAME_2]
        res_num2 = row[RES_NUM_2]
        chain2 = row[CHAIN_2]
        atom2 = row[ATOM_2]
        atom_1_coords = [row["x1"], row["y1"], row["z1"]]
        atom_2_coords = [row["x2"], row["y2"], row["z2"]]
        rings_to_check_for_atom1 = df.loc[df[RES_NUM_1] == res_num1]
        rings_to_check_for_atom2 = df.loc[df[RES_NUM_2] == res_num2]
        for temp_index1, temp_row1 in rings_to_check_for_atom1.iterrows():
            if temp_row1[RES_NAME_2] in residues:
                if res_name1 in lone_pair_pi and atom1 in lone_pair_pi[res_name1]:
                    check_lone_pair_pi(atom_1_coords, row, temp_row1, ring_params, True)
                if res_name1 in anion_pi and atom1 in anion_pi[res_name1]:
                    check_anion_pi(atom_1_coords, row, temp_row1, ring_params, True)
                if res_name1 in cation_pi and atom1 in cation_pi[res_name1]:
                    check_cation_pi(atom_1_coords, row, temp_row1, ring_params, True)
                if res_name1 in amino_pi and atom1 in amino_pi[res_name1]:
                    check_amino_pi(atom_1_coords, row, temp_row1, ring_params, True)
        for temp_index2, temp_row2 in rings_to_check_for_atom2.iterrows():
            if temp_row2[RES_NAME_1] in residues:
                if res_name2 in lone_pair_pi and atom2 in lone_pair_pi[res_name2]:
                    check_lone_pair_pi(
                        atom_2_coords, row, temp_row2, ring_params, False
                    )
                if res_name2 in anion_pi and atom2 in anion_pi[res_name2]:
                    check_anion_pi(atom_2_coords, row, temp_row2, ring_params, False)
                if res_name2 in cation_pi and atom2 in cation_pi[res_name2]:
                    check_cation_pi(atom_2_coords, row, temp_row2, ring_params, False)
                if res_name2 in amino_pi and atom2 in amino_pi[res_name2]:
                    check_amino_pi(atom_2_coords, row, temp_row2, ring_params, False)
    return (
        lone_pair_pi_interactions,
        anion_pi_interactions,
        cation_pi_interactions,
        amino_pi_interactions,
    )
