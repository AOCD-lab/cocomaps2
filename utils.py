"""
Comprehensive utility module for processing and analyzing biomolecular structures.
Includes functions for parsing PDB/mmCIF files, computing atomic distances, angles,
and planes, and identifying geometric or interaction-based features (e.g., π-systems,
metal-mediated contacts, hydrogen bonds, and van der Waals interactions).  

Provides helper methods for building per-chain and per-residue datasets, filtering
atoms, merging connected rings, generating printable PDB lines, and constructing
distance tables between chains.  

Also includes tools for reformatting and exporting PDB files — such as complex,
chain-specific, and truncated versions — while maintaining proper atom alignment
and chain mapping for downstream computational analyses.
"""

import numpy as np
import pandas as pd
import math
from scipy.spatial import ConvexHull, Delaunay
from bidict import bidict
import numpy as np
import copy
from Bio.PDB import PDBParser, PDBIO, Select

from update_constants import tc
from constants import METAL_CUTOFF

CUT_OFF = tc.CUT_OFF
APOLAR_DIST = tc.APOLAR_TOLERANCE
from constants import (
    CHAIN1_IDENTIFIER,
    CHAIN2_IDENTIFIER,
    pdb_fromat_dict,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_NAME,
    WATER_MEDIATED_NAME,
    WATER_IDENTITY_NAME,
    METAL_MEDIATED_NAME,
    METAL_IDENTITY_NAME,
    ATOM_1_ORG_NAME,
    ATOM_2_ORG_NAME,
    HYDROGEN_BOND_DICT_STRUCTURE,
    SALT_BRIDGE_DICT_STRUCTURE
)
from read_pdb_regex import read_pdb
from collections import defaultdict
from residues import (
    residues,
    C_H_PI,
    NSOH_PI,
    METAL_MEDIATED,
    METALS,
    WEAK_H_BONDS,
    ACCEPTOR_ATOMS,
    RADII,
)


def def_value1():
    return {}


def def_value2():
    return {"C_line": None, "H_lines": list}


def reverse_mapping_chains_in_df(chain_mapping, df, interaction_type=None):
    if interaction_type == WATER_MEDIATED_NAME:
        for index, row in df.iterrows():
            try:
                water_id_chain = row[WATER_IDENTITY_NAME].split("_")[0]
                if water_id_chain in chain_mapping:
                    water_id_chain = chain_mapping[water_id_chain]
                    df.loc[index, WATER_IDENTITY_NAME] = (
                        f"{water_id_chain}_{row[WATER_IDENTITY_NAME].split('_')[1]}_HOH"
                    )
                df.loc[index, CHAIN_1] = chain_mapping[row[CHAIN_1]]
                df.loc[index, CHAIN_2] = chain_mapping[row[CHAIN_2]]
            except Exception as e:
                print(e)
        return df
    elif interaction_type == METAL_MEDIATED_NAME:
        for index, row in df.iterrows():
            try:
                metal_id_chain = row[METAL_IDENTITY_NAME].split("_")[0]
                if metal_id_chain in chain_mapping:
                    metal_id_chain = chain_mapping[metal_id_chain]
                    df.loc[index, METAL_IDENTITY_NAME] = (
                        f"{metal_id_chain}_{row[METAL_IDENTITY_NAME].split('_')[1]}_HOH"
                    )
                df.loc[index, CHAIN_1] = chain_mapping[row[CHAIN_1]]
                df.loc[index, CHAIN_2] = chain_mapping[row[CHAIN_2]]
            except Exception as e:
                print(e)
        return df
    else:
        for index, row in df.iterrows():
            try:
                df.loc[index, CHAIN_1] = chain_mapping[row[CHAIN_1]]
                df.loc[index, CHAIN_2] = chain_mapping[row[CHAIN_2]]
            except Exception as e:
                print(e)
        return df


def get_combined_rings(all_rings):
    ind_dict = {}
    iterator = 0

    for i in range(len(all_rings)):
        ind_dict[i] = []
    for ind1, ring1 in enumerate(all_rings):
        iterator += 1
        for ind2, ring2 in enumerate(all_rings[iterator:]):
            ind2 += iterator
            if len(set(ring1).intersection(ring2)) >= 2:
                ind_dict[ind1].append(ind2)
    ind_dict

    connections = []
    for key, value in ind_dict.items():
        if not connections:
            if value:
                connections.append([val for val in value])
                connections[0].append(key)
        else:
            if value:
                values = [val for val in value]
                values.append(key)
                connections.append(values)

    full_rings = []
    added_inds = []

    for connection in connections:
        temp_ring = []
        for ring_ind in connection:
            added_inds.append(ring_ind)
            temp_ring += all_rings[ring_ind]
        full_rings.append(temp_ring)

    all_inds = set([i for i in range(len(all_rings))])
    added_inds = set(added_inds)
    remaining_inds = all_inds.difference(added_inds)

    for remaining_ind in remaining_inds:
        full_rings.append(all_rings[remaining_ind])
    for ind, ring in enumerate(full_rings):
        full_rings[ind] = list(set(ring))

    return full_rings


def get_distance_table(focus_lines):
    threshold = CUT_OFF
    A = focus_lines[list(focus_lines.keys())[0]]  # Your list A
    B = focus_lines[list(focus_lines.keys())[1]]  # Your list B
    metal_oriented_lines = {
        list(focus_lines.keys())[0]: [],
        list(focus_lines.keys())[1]: [],
    }
    # Extract coordinates and other necessary information
    coords1 = np.array(
        [[float(x), float(y), float(z)] for x, y, z in (residue[6:9] for residue in A)]
    )
    coords2 = np.array(
        [[float(x), float(y), float(z)] for x, y, z in (residue[6:9] for residue in B)]
    )

    # Compute pairwise distances using broadcasting
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=-1))

    # Mask for distances within the threshold
    mask = distances <= threshold
    mask2 = distances <= METAL_CUTOFF

    distance_table = []
    for i, residue1 in enumerate(A):
        for j, residue2 in enumerate(B):
            if mask[i, j]:
                dist = round(distances[i, j], 2)
                line = [
                    residue1[3],
                    residue1[5],
                    residue1[4],
                    residue1[2],
                    coords1[i][0],
                    coords1[i][1],
                    coords1[i][2],
                    residue2[3],
                    residue2[5],
                    residue2[4],
                    residue2[2],
                    coords2[j][0],
                    coords2[j][1],
                    coords2[j][2],
                    dist,
                    residue1[11],
                    residue2[11],
                ]
                distance_table.append(line)
            if mask2[i, j]:
                metal_oriented_lines[list(focus_lines.keys())[0]].append(residue1)
                metal_oriented_lines[list(focus_lines.keys())[1]].append(residue2)
    columns = [
        RES_NAME_1,
        RES_NUM_1,
        CHAIN_1,
        ATOM_1,
        "x1",
        "y1",
        "z1",
        RES_NAME_2,
        RES_NUM_2,
        CHAIN_2,
        ATOM_2,
        "x2",
        "y2",
        "z2",
        DISTANCE_NAME,
        ATOM_1_ORG_NAME,
        ATOM_2_ORG_NAME
        
    ]
    df = pd.DataFrame(distance_table)
    df.columns = columns
    full_df = df.copy()
    full_df.to_csv("full.csv")

    return full_df, metal_oriented_lines


def get_distance_table_(focus_lines):
    threshold = CUT_OFF
    A = focus_lines[list(focus_lines.keys())[0]]  # Your list A
    B = focus_lines[list(focus_lines.keys())[1]]  # Your list B

    # Preallocate distance_table
    distance_table = []

    for i, residue1 in enumerate(A):
        x1, y1, z1 = map(float, residue1[6:9])

        for j, residue2 in enumerate(B):
            x2, y2, z2 = map(float, residue2[6:9])

            dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
            dist = round(dist, 2)

            if dist <= threshold:
                line = [
                    residue1[3],
                    residue1[5],
                    residue1[4],
                    residue1[2],
                    x1,
                    y1,
                    z1,
                    residue2[3],
                    residue2[5],
                    residue2[4],
                    residue2[2],
                    x2,
                    y2,
                    z2,
                    dist,
                ]
                distance_table.append(line)

    columns = [
        RES_NAME_1,
        RES_NUM_1,
        CHAIN_1,
        ATOM_1,
        "x1",
        "y1",
        "z1",
        RES_NAME_2,
        RES_NUM_2,
        CHAIN_2,
        ATOM_2,
        "x2",
        "y2",
        "z2",
        DISTANCE_NAME,
    ]
    df = pd.DataFrame(distance_table)
    df.columns = columns
    full_df = df.copy()
    full_df.to_csv("full.csv")

    return full_df


def get_constituents_from_line(parsed_line):
    atom_name = parsed_line[2]
    res_name = parsed_line[3]
    res_num = parsed_line[5]
    chain = parsed_line[4]
    x = parsed_line[6]
    y = parsed_line[7]
    z = parsed_line[8]
    return atom_name, res_name, res_num, chain, x, y, z


def get_residue_chunks(full_pdb_lines, chains):
    chunks = {}
    for chain in chains:
        chunks[chain] = {}
    for line in full_pdb_lines:
        if line[4] in chains:
            if line[5] not in chunks[line[4]].keys():
                chunks[line[4]][line[5]] = []
            chunks[line[4]][line[5]].append(line)
    return chunks


def get_closest_atom_in_residue(focus_atom_line, residue_chunk):
    closest_atom_line = None
    dist = 100
    for line in residue_chunk:
        check_dist = get_distance_between_lines(focus_atom_line, line)
        if check_dist < dist and line[2] != focus_atom_line[2]:
            closest_atom_line = line
            dist = check_dist
    return closest_atom_line


def get_items_from_identity(identity):
    ind = None
    for i, letter in enumerate(identity[::-1]):
        if not letter.isdigit():
            ind = i
            break
    if ind:
        b = identity[:-ind]
        c = identity[-ind:]
        return b, c
    else:
        return None, None


def get_c_h_n_h_lines(chains_lines):
    c_h_combos = {}
    n_h_combos = {}
    for chain, chain_lines in chains_lines.items():
        for chain_line in chain_lines:
            c_h_combos = add_c_h_lines_to_dict(
                c_h_combos, chain_line, chains_lines.keys()
            )
            n_h_combos = add_n_h_lines_to_dict(
                n_h_combos, chain_line, chains_lines.keys()
            )
    return c_h_combos, n_h_combos


def add_c_h_lines_to_dict(c_h_combos, line, rna_chains):
    atom_name, res_name, res_num, chain, x, y, z = get_constituents_from_line(line)
    if chain in rna_chains:
        if chain not in c_h_combos:
            c_h_combos[chain] = {}
        if res_name in C_H_PI and atom_name in C_H_PI[res_name]:
            if res_num not in c_h_combos[chain]:
                c_h_combos[chain][res_num] = {}
            if atom_name not in c_h_combos[chain][res_num]:
                c_h_combos[chain][res_num][atom_name] = {
                    "c_coords": [x, y, z],
                    "h_atoms": [],
                }
        elif res_name in C_H_PI and atom_name.startswith("H"):
            for c_atom, h_atoms in C_H_PI[res_name].items():
                if (
                    res_num in c_h_combos[chain]
                    and c_atom in c_h_combos[chain][res_num]
                    and atom_name in h_atoms
                ):
                    c_h_combos[chain][res_num][c_atom]["h_atoms"].append(
                        {"name": atom_name, "coords": [x, y, z]}
                    )
    return c_h_combos


def add_n_h_lines_to_dict(n_h_combos, line, rna_chains):
    atom_name, res_name, res_num, chain, x, y, z = get_constituents_from_line(line)
    if chain in rna_chains:
        if chain not in n_h_combos:
            n_h_combos[chain] = {}
        if res_name in NSOH_PI and atom_name in NSOH_PI[res_name]:
            if res_num not in n_h_combos[chain]:
                n_h_combos[chain][res_num] = {}
            if atom_name not in n_h_combos[chain][res_num]:
                n_h_combos[chain][res_num][atom_name] = {
                    "n_coords": [x, y, z],
                    "h_atoms": [],
                }
        elif res_name in NSOH_PI and atom_name.startswith("H"):
            for n_atom, h_atoms in NSOH_PI[res_name].items():
                if (
                    res_num in n_h_combos
                    and n_atom in n_h_combos[chain][res_num]
                    and atom_name in h_atoms
                ):
                    n_h_combos[chain][res_num][n_atom]["h_atoms"].append(
                        {"name": atom_name, "coords": [x, y, z]}
                    )
    return n_h_combos


def get_metal_and_acceptor_lines(
    all_parsed_lines, chains_mapping: bidict, metal_oriented_lines
):
    """This function in general provides with lines from the pdb file which
    include either a metal or a acceptor atom, acceptor atoms are determined by METAL_MEDIATED
    dictionary in ressidues file.
    param pdb_file: name or path of the pdb file
    param chains: list of chains. Example ['A','B']
    Return metal_atoms_lines, donor_atom_lines: lists of such lines"""

    metal_atoms_lines = []
    acceptor_atom_lines_rna = {CHAIN1_IDENTIFIER: [], CHAIN2_IDENTIFIER: []}

    for line in all_parsed_lines:
        line[11] = ''.join(filter(str.isalpha, line[11]))
        if line[0] == "HETATM" and ''.join(filter(str.isalpha, line[11])) in METALS:
            if line[4] in chains_mapping.keys():
                line[4] = chains_mapping[line[4]]
            metal_atoms_lines.append(line)
    for chain, chain_lines in metal_oriented_lines.items():
        for line in chain_lines:
            if line[3] in METAL_MEDIATED.keys() and line[2] in METAL_MEDIATED[line[3]]:
                acceptor_atom_lines_rna[chain].append(line)

    return metal_atoms_lines, acceptor_atom_lines_rna


def make_rings(chain_lines):
    rings = {}
    for chain, parsed_lines in chain_lines.items():
        for line in parsed_lines:
            atom_name, res_name, res_num, chain_, x, y, z = get_constituents_from_line(
                line
            )
            if chain not in rings:
                rings[chain] = {}
            if res_name in residues:
                atoms = residues[res_name]
                if res_num in rings[chain] and atom_name in atoms:
                    rings[chain][res_num][atom_name] = [x, y, z]
                elif res_num not in rings[chain] and atom_name in atoms:
                    res_atom_dict = {}
                    for i in atoms:
                        res_atom_dict[i] = []
                    rings[chain][res_num] = res_atom_dict
                    rings[chain][res_num][atom_name] = [x, y, z]
    return rings


def initial_filter_data_frame_creator(
    chains_lines,
):

    # distance_table = distance_table_creator(chains)
    rings = make_rings(chains_lines)
    df, metal_oriented_lines = get_distance_table(chains_lines)
    # columns = [
    #     RES_NAME_1,
    #     RES_NUM_1,
    #     CHAIN_1,
    #     ATOM_1,
    #     "x1",
    #     "y1",
    #     "z1",
    #     RES_NAME_2,
    #     RES_NUM_2,
    #     CHAIN_2,
    #     ATOM_2,
    #     "x2",
    #     "y2",
    #     "z2",
    #     DISTANCE_NAME,
    # ]
    # df = pd.DataFrame(distance_table)
    # df.columns = columns
    full_df = df.copy()
    full_df.to_csv("trial24.csv")
    full_df = full_df.drop_duplicates()
    df = df[~df[ATOM_1].str.startswith("H")]
    df = df[~df[ATOM_2].str.startswith("H")]
    df = (
        df.sort_values("Distance (Å)", ascending=True)
        .drop_duplicates(
            subset=[
                RES_NAME_1,
                RES_NUM_1,
                CHAIN_1,
                RES_NAME_2,
                RES_NUM_2,
                CHAIN_2,
            ]
        )
        .sort_index()
    )
    filetered_df = df[df[RES_NAME_1] != "HOH"]
    filetered_df = filetered_df[filetered_df[RES_NAME_2] != "HOH"]
    # df.to_csv("trial23.csv")
    return full_df, filetered_df, rings, metal_oriented_lines


def get_atom_name_res_num_from_index_ligand(ligand_parsed_lines, index):
    return (
        ligand_parsed_lines[index][2],
        ligand_parsed_lines[index][3],
        ligand_parsed_lines[index][5],
    )


def get_angle_between_3_lines(line1, line2, line3):
    a = np.array([line1[6], line1[7], line1[8]])
    b = np.array([line2[6], line2[7], line2[8]])
    c = np.array([line3[6], line3[7], line3[8]])
    return get_angle_between_3_points(a, b, c)


def get_angle_between_3_points(p1, p2, p3):
    a = np.array(p1)
    b = np.array(p2)
    c = np.array(p3)

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


def check_apolar_vdw(full_df: pd.DataFrame):
    for index, row in full_df.iterrows():
        try:
            # self.check_s_bonds(row)
            # # self.check_vdw_bonds(row)
            # self.check_clash_bonds(row)
            res_1 = row["RNA Res. Name"]
            res_2 = row["Ligand Res. Name"]
            hydrobhobic_atom1 = row["RNA Atom"]
            hydrobhobic_atom2 = row["Ligand Atom"]
            if (
                hydrobhobic_atom2 == "C5"
                and res_2 == "U"
                and row["Distance (Å)"] <= APOLAR_DIST
            ):
                pass  # add to main table
            elif (
                hydrobhobic_atom2 == "C5"
                and res_2 == "C"
                and row["Distance (Å)"] <= APOLAR_DIST
            ):
                pass  # add to main table
        except:
            pass


def check_similar_rows_between_dicsts(dict1, dict2):

    # Find matching rows across all possible indices
    matching_rows = []

    for i in range(len(dict1[list(dict1.keys())[0]])):
        for j in range(len(dict2[list(dict1.keys())[0]])):
            if all(str(dict1[key][i]) == str(dict2[key][j]) for key in dict1.keys()):
                matching_rows.append((i, j))

    return matching_rows


def remove_similar_rows_from_dicts(dict_to_check, dict_from_check, dict_structure):
    matching_rows = check_similar_rows_between_dicsts(dict_to_check, dict_from_check)
    filtered_dict = copy.deepcopy(dict_structure)
    indexes_to_remove = [i for i, j in matching_rows]
    for key, value in dict_to_check.items():
        for ind, val in enumerate(value):
            if ind not in indexes_to_remove:
                filtered_dict[key].append(val)
    return filtered_dict

def filter_hydrogen_salt_dicts(dict_hydrogen, dict_salt):
    hydrogen_indexs_to_remove = []
    salt_ind, hydrogen_ind = 0,0
    for res1_salt, res2_salt in zip(dict_salt[RES_NUM_1], dict_salt[RES_NUM_2]):
        hydrogen_ind = 0
        for res1_hydrogen, res2_hydrogen in zip(dict_hydrogen[RES_NUM_1], dict_hydrogen[RES_NUM_2]):
            if res1_hydrogen == res1_salt and res2_hydrogen == res2_salt:
                hydrogen_indexs_to_remove.append(hydrogen_ind)
            hydrogen_ind +=1
        salt_ind+=1
    n = len(dict_hydrogen[RES_NUM_1])
    indices_to_keep = [i for i in range(n) if i not in hydrogen_indexs_to_remove]
    for key in dict_hydrogen:
        dict_hydrogen[key] = [dict_hydrogen[key][i] for i in indices_to_keep]
    
    return dict_hydrogen

def remove_salt_duplicates(salt_bridge_dict):
    indices_to_remove = set()
    n = len(salt_bridge_dict[RES_NUM_1])
    
    # Iterate over all pairs of indices
    for i in range(n):
        for j in range(i + 1, n):
            if (
                salt_bridge_dict[RES_NUM_1][i] == salt_bridge_dict[RES_NUM_1][j] and
                salt_bridge_dict[RES_NUM_2][i] == salt_bridge_dict[RES_NUM_2][j]
            ):
                # Compare DISTANCE_NAME values and mark the larger one for removal
                if salt_bridge_dict[DISTANCE_NAME][i] > salt_bridge_dict[DISTANCE_NAME][j]:
                    indices_to_remove.add(i)
                else:
                    indices_to_remove.add(j)
    
    # Remove marked indices
    indices_to_keep = [i for i in range(n) if i not in indices_to_remove]
    for key in salt_bridge_dict:
        salt_bridge_dict[key] = [salt_bridge_dict[key][i] for i in indices_to_keep]
    
    return salt_bridge_dict
    




def get_heavy_atoms_from_groups(groups):
    combos = {}
    for group in groups:
        heavy_atom = int(group.split("-")[0])
        h_atom = int(group.split("-")[-1])
        if heavy_atom not in combos:
            combos[heavy_atom] = []
        combos[heavy_atom].append(h_atom)
    return combos


def create_printable_line(line):
    printable_line = []
    for i in range(79):
        printable_line.append(" ")
    for ind, element in enumerate(line):
        element = list(str(element))
        fomat = pdb_fromat_dict[ind]
        if fomat[1] == "l":
            if ind == 1 and len(element) > 4:
                element = list(element)
                element = element[:4]
                element = "".join(element)
            printable_line[fomat[0][0] : len(element) + fomat[0][0] + 1] = element
        else:
            if ind == 1 and len(element) > 4:
                element = list(element)
                element = element[:5]
            printable_line[fomat[0][1] - len(element) + 1 : fomat[0][1] + 1] = element
    return "".join(printable_line)


def create_plane(point1, point2, point3):
    # Find two vectors lying on the plane
    vector1 = np.array(point2) - np.array(point1)
    vector2 = np.array(point3) - np.array(point1)

    # Find the normal vector of the plane
    normal_vector = np.cross(vector1, vector2)

    # Use the normal vector and one of the points to form the equation of the plane
    a, b, c = normal_vector
    deviding_factor = (math.sqrt(a**2+b**2+c**2))
    a1 = a/deviding_factor
    b1 = b/deviding_factor
    c1 = c/deviding_factor
    d = -(a * point1[0] + b * point1[1] + c * point1[2])
    d1 = d/deviding_factor
    # Return the coefficients of the plane equation
    return a1, b1, c1, d1


def project_point_onto_plane(point, plane_coefficients):
    a, b, c, d = plane_coefficients

    # Calculate the distance from the point to the plane
    k = (-d - a * point[0] - b * point[1] - c * point[2]) #/ ((a**2 + b**2 + c**2))
    xp = point[0] + k * a
    yp = point[1] + k * b
    zp = point[2] + k * c
    # Calculate the projected point

    return (xp, yp, zp), np.abs(k) #this is try


def get_centroid(points):
    x, y, z = 0, 0, 0
    for point in points:
        x += point[0]
        y += point[1]
        z += point[2]
    return [x / len(points), y / len(points), z / len(points)]


def get_dist_between_points(p1, p2):
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    z = p1[2] - p2[2]
    d = math.sqrt(x * x + y * y + z * z)
    return d


def get_distance_between_lines(line1, line2):
    """This function gets the distance between 2 pdb ATOM lines
    param line1, line3: list of a pdb line with one part of line as one element. Example ['ATOM','C','OP1'......]
    Return dist: float distance thus calculated"""

    x1 = float(line1[6])
    y1 = float(line1[7])
    z1 = float(line1[8])
    x2 = float(line2[6])
    y2 = float(line2[7])
    z2 = float(line2[8])
    a = np.array((x1, y1, z1))
    b = np.array((x2, y2, z2))
    dist = np.linalg.norm(a - b)

    return dist


def is_point_inside_3d_polygon(point, polygon):
    hull = ConvexHull(polygon)
    new_hull = ConvexHull(np.concatenate((polygon, [point])))
    return np.array_equal(new_hull.vertices, hull.vertices)


def get_angle_atom_h_centroid(c_coord, h_coord, centroid_coord):
    a = np.array(c_coord)
    b = np.array(h_coord)
    c = np.array(centroid_coord)

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


def get_c_centoid_dist(c_coords, centroid):
    return get_dist_between_points(c_coords, centroid)


def get_angle_normal_atom(c_coords, centroid_coord, plane_coefficients, r_cen):
    centroid_vector = np.array(c_coords) - np.array(centroid_coord)
    a1, b1, c1, d1 = plane_coefficients
    normal_to_ring = np.array([a1, b1, c1])
    dot_prod = np.linalg.norm(np.dot(normal_to_ring, centroid_vector))
    cos_theta = dot_prod / (r_cen * math.sqrt(a1 * a1 + b1 * b1 + c1 * c1))
    return np.degrees(np.arccos(cos_theta))


def check_atom_h_pi_condition(
    c_coord,
    h_coord,
    plane_coefficients,
    centroid,
    dist_threshold,
    theta1_threshold,
    theta2_threshold,
):
    r_cen = get_c_centoid_dist(c_coord, centroid)
    atom_h_centoid_angle = get_angle_atom_h_centroid(c_coord, h_coord, centroid)
    atom_normal_angle = get_angle_normal_atom(
        c_coord, centroid, plane_coefficients, r_cen
    )

    if (
        r_cen <= dist_threshold
        and atom_h_centoid_angle >= theta1_threshold
        and atom_normal_angle <= theta2_threshold
    ):
        return r_cen, atom_h_centoid_angle, atom_normal_angle
    else:
        return None, None, None


def get_constituents_from_row(row):
    res_name1 = row[RES_NAME_1]
    atom1 = row[ATOM_1]
    res_name2 = row[RES_NAME_2]
    atom2 = row[ATOM_2]
    res_num1 = row[RES_NUM_1]
    res_num2 = row[RES_NUM_2]
    chain1 = row[CHAIN_1]
    chain2 = row[CHAIN_2]
    x1, y1, z1 = row["x1"], row["y1"], row["z1"]
    x2, y2, z2 = row["x2"], row["y2"], row["z2"]
    return (
        res_name1,
        atom1,
        res_name2,
        atom2,
        res_num1,
        res_num2,
        chain1,
        chain2,
        [x1, y1, z1],
        [x2, y2, z2],
    )


def check_clash_in_row(row):
    try:
        atom1 = row[ATOM_1]
        atom2 = row[ATOM_2]
        dist = row[DISTANCE_NAME]
        radii_dist = RADII[atom1[0]] + RADII[atom2[0]]
        if float(dist) < radii_dist:
            return True
        return False
    except Exception as e:
        print(e)
        return False


def get_acceptor_atom_and_hbond_lines(chains_parsed_lines, chains):
    weak_h_bond_lines = {}
    for chain in chains:
        weak_h_bond_lines[chain] = defaultdict(def_value1)
    acceptor_atom_lines = {chains[0]: [], chains[1]: []}
    current_res_num = chains_parsed_lines[0][5]
    for line in chains_parsed_lines:
        if line[3] in WEAK_H_BONDS.keys() and line[4] in chains:
            temp_res_num = line[5]
            temp_res_name = line[3]
            chain = line[4]
            if (
                line[2] in WEAK_H_BONDS[line[3]].keys()
                and temp_res_num == current_res_num
            ):
                weak_h_bond_lines[chain][temp_res_num][line[2]] = defaultdict(
                    def_value2
                )
                weak_h_bond_lines[chain][temp_res_num][line[2]]["C_line"] = line
                weak_h_bond_lines[chain][temp_res_num][line[2]]["H_lines"] = []
            for c_atom, h_atoms in WEAK_H_BONDS[line[3]].items():
                if (
                    line[2] in h_atoms
                    and c_atom in weak_h_bond_lines[chain][temp_res_num]
                ):
                    weak_h_bond_lines[chain][temp_res_num][c_atom]["H_lines"].append(
                        line
                    )
            current_res_num = temp_res_num

        if (
            line[3] in ACCEPTOR_ATOMS.keys()
            and line[2] in ACCEPTOR_ATOMS[line[3]]
            and line[4] in chains
        ):
            acceptor_atom_lines[line[4]].append(line)

    return weak_h_bond_lines, acceptor_atom_lines


def create_printable_line_cif_version(line):
    printable_line = []
    for i in range(79):
        printable_line.append(" ")
    for ind, element in enumerate(line):
        element = list(str(element))
        fomat = pdb_fromat_dict[ind]
        if fomat[1] == "l":
            if ind == 1 and len(element) > 4:
                element = list(element)
                element = element[:4]
                element = "".join(element)
            if ind in [6,7,8]: #####
                formatted = f"{element:.3f}"
                element = formatted
            printable_line[fomat[0][0] : len(element) + fomat[0][0] + 1] = element
        else:
            if ind == 1 and len(element) > 4:
                element = list(element)
                element = element[:5]
            printable_line[fomat[0][1] - len(element) + 1 : fomat[0][1] + 1] = element
    return "".join(printable_line)


def get_chains_ranges_lines(
    parsed_lines: list,
    ranges1: list,
    ranges2: list,
    chains,
    pdb_path: str,
):
    chains_lines = {CHAIN1_IDENTIFIER: [], CHAIN2_IDENTIFIER: []}
    chain1 = chains[0]
    chain2 = chains[1]
    chains_mappings = {CHAIN1_IDENTIFIER: chain1, CHAIN2_IDENTIFIER: chain2}
    # chains_mappings = bidict(chains_mappings)
    all_pdb_changed_chains_file_path = f"{pdb_path}_formatted.pdb"
    all_pdb_changed_chains_file = open(all_pdb_changed_chains_file_path, "w")
    compound_pdb_path = f"{pdb_path}_{chain1}_{chain2}_complex.pdb"
    chain1_pdb_path = f"{pdb_path}_{chain1}.pdb"
    chain2_pdb_path = f"{pdb_path}_{chain2}.pdb"
    compound_pdb = open(compound_pdb_path, "w")
    chain1_pdb = open(chain1_pdb_path, "w")
    chain2_pdb = open(chain2_pdb_path, "w")
    all_pdb_parsed_lines = []
    all_pdb_raw_lines = []
    chain1_line_written = False
    chain2_line_written = False
    for line in parsed_lines:
        line_chain = line[4]
        line_res_num = line[5]
        if line_chain == chain1:
            new_line1 = copy.deepcopy(line)
            if (
                ranges1
                and line_res_num.isdigit()
                and int(line_res_num) >= ranges1[0]
                and int(line_res_num) <= ranges1[1]
            ):
                new_line1[4] = CHAIN1_IDENTIFIER
                chains_lines[CHAIN1_IDENTIFIER].append(new_line1)
                printable_line = create_printable_line(new_line1)
                compound_pdb.write(f"{printable_line}\n")
                chain1_pdb.write(f"{printable_line}\n")
            elif not ranges1:

                new_line1[4] = CHAIN1_IDENTIFIER
                chains_lines[CHAIN1_IDENTIFIER].append(new_line1)
                printable_line = create_printable_line(new_line1)
                compound_pdb.write(f"{printable_line}\n")
                chain1_pdb.write(f"{printable_line}\n")
            else:
                printable_line = create_printable_line(new_line1)
            all_pdb_parsed_lines.append(new_line1)
            all_pdb_changed_chains_file.write(f"{printable_line}\n")
            chain1_line_written = True

        if line_chain == chain2:
            new_line2 = copy.deepcopy(line)
            if (
                ranges2
                and line_res_num.isdigit()
                and int(line_res_num) >= ranges2[0]
                and int(line_res_num) <= ranges2[1]
            ):
                new_line2[4] = CHAIN2_IDENTIFIER
                chains_lines[CHAIN2_IDENTIFIER].append(new_line2)
                printable_line = create_printable_line(new_line2)
                compound_pdb.write(f"{printable_line}\n")
                chain2_pdb.write(f"{printable_line}\n")
            elif not ranges2:
                new_line2[4] = CHAIN2_IDENTIFIER
                chains_lines[CHAIN2_IDENTIFIER].append(new_line2)
                printable_line = create_printable_line(new_line2)
                compound_pdb.write(f"{printable_line}\n")
                chain2_pdb.write(f"{printable_line}\n")
            else:
                printable_line = create_printable_line(new_line2)
            all_pdb_parsed_lines.append(new_line2)
            all_pdb_changed_chains_file.write(f"{printable_line}\n")
            chain2_line_written = True
        if not chain1_line_written and not chain2_line_written:
            all_pdb_parsed_lines.append(line)
            printable_line = create_printable_line(line)
            all_pdb_changed_chains_file.write(f"{printable_line}\n")
        chain1_line_written, chain2_line_written = False, False
    all_pdb_changed_chains_file.close()
    compound_pdb.close()
    chain1_pdb.close()
    chain2_pdb.close()
    return (
        chains_lines,
        chains_mappings,
        all_pdb_parsed_lines,
        all_pdb_changed_chains_file_path,
        compound_pdb_path,
        chain1_pdb_path,
        chain2_pdb_path,
    )

def chain_selection_pdbs_making(chain1, chain2, pdb_path):

    # Load the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("molecule", pdb_path)
    compound_pdb_path = f"{pdb_path}_{chain1}_{chain2}_complex.pdb"
    chain1_pdb_path = f"{pdb_path}_{chain1}.pdb"
    chain2_pdb_path = f"{pdb_path}_{chain2}.pdb"

    io = PDBIO()

    # Define a chain selector
    class ChainSelect(Select):
        def __init__(self, chains):
            self.chains = set(chains)
        def accept_chain(self, chain):
            return chain.get_id() in self.chains

    # Save only chain A
    io.set_structure(structure)
    io.save(chain1_pdb_path, ChainSelect([chain1]))

    # Save only chain B
    io.set_structure(structure)
    io.save(chain2_pdb_path, ChainSelect([chain2]))

    # Save chains A and B together
    io.set_structure(structure)
    io.save(compound_pdb_path, ChainSelect([chain1, chain2]))
    return chain1_pdb_path, chain2_pdb_path, compound_pdb_path


def make_truncated_pdb(water_idetities, metal_identities, chains_mapping, all_parsed_lines, pdb_path):
    truncated_pdb_path = pdb_path.replace("complex", "truncated")
    truncated_pdb_file = open(truncated_pdb_path, "w")
    water_ids = [(i.split("_")[0], str(int(i.split("_")[1])), i.split("_")[2]) for i in water_idetities if len(i.split("_")) == 3]
    metal_ids = [(i.split("_")[0], str(int(i.split("_")[1])), i.split("_")[2]) for i in metal_identities if len(i.split("_")) == 3]
    for line in all_parsed_lines:
        line_id = (line[4], line[5], line[3])
        if line[4] in chains_mapping and line[3] != "HOH":
            new_line = copy.deepcopy(line)
            new_line[4] = chains_mapping[line[4]]
            printable_line = create_printable_line_cif_version(new_line)
            truncated_pdb_file.write(f"{printable_line}\n")
        elif line_id in water_ids:
            new_line = copy.deepcopy(line)
            if line_id[0] in chains_mapping:
                new_line[4] = chains_mapping[line_id[0]]
            printable_line = create_printable_line_cif_version(new_line)
            truncated_pdb_file.write(f"{printable_line}\n")
        elif line_id in metal_ids:
            new_line = copy.deepcopy(line)
            if line_id[0] in chains_mapping:
                new_line[4] = chains_mapping[line_id[0]]
            printable_line = create_printable_line_cif_version(new_line)
            truncated_pdb_file.write(f"{printable_line}\n")
    truncated_pdb_file.close()
        



