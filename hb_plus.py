"""
Runs HBPLUS and HBADD to detect hydrogen bonds and water-mediated interactions
from PDB structures. Parses output files, classifies direct and water-bridged
hydrogen bonds between specified chains, and returns structured dictionaries
of interaction data for further analysis.
"""

import os
from collections import defaultdict
import copy
from constants import (
    HB_PLUS_PATH,
    HB_PLUS_FISRT_LINE,
    HYDROGEN_BOND_DICT_STRUCTURE,
    WATER_MEDIATED_DICT_STRUCTURE,
    HBADD_PATH,
    HBADD_EMPTY_FILE_PATH,
    HBADD_FILE_NAME,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_NAME,
    WATER_IDENTITY_NAME,
    DISTANCE_FROM_RES1_NAME,
    DISTANCE_FROM_RES2_NAME,
    DHA_ANGLE_NAME,
)
from update_constants import tc
from utils import get_items_from_identity, remove_similar_rows_from_dicts

HBOND_DIST = tc.HBOND_DIST
HBOND_ANGLE = tc.HBOND_ANGLE
SBRIDGE_DIST = tc.SBRIDGE_DIST
WBRIDGE_DIST = tc.WBRIDGE_DIST

HB_PLUS_DIST = max([HBOND_DIST, WBRIDGE_DIST, SBRIDGE_DIST])


def def_value():
    return []


def run_hbadd(pdb_file_path):
    head, tail = os.path.split(pdb_file_path)
    os.system(f"{HBADD_PATH} {pdb_file_path} {HBADD_EMPTY_FILE_PATH} -wkdir {head}")
    hbadd_file = os.path.join(head, HBADD_FILE_NAME)
    return hbadd_file


def run_hbplus(pdb_file_path):
    hbadd_file = run_hbadd(pdb_file_path)
    hbplus_command = f"{HB_PLUS_PATH} -D {HB_PLUS_DIST} -a {HBOND_ANGLE} -f {hbadd_file} {pdb_file_path}"
    os.system(hbplus_command)
    # os.system(f"{HB_PLUS_PATH} -D {HB_PLUS_DIST} -a {HBOND_ANGLE}  {pdb_file_path}")
    head, tail = os.path.split(pdb_file_path)
    current_file_name = f'{".".join(tail.split(".")[:-1])}.hb2'
    hbplus_file = os.path.join(os.getcwd(), current_file_name)
    if not os.path.isfile(hbplus_file):
        os.system(f"{HB_PLUS_PATH} -D {HB_PLUS_DIST} -a {HBOND_ANGLE} {pdb_file_path}")
    return hbplus_file


def get_hoh_and_h_lines(line_constituents):
    HOH_lines = defaultdict(def_value)
    H_lines = {}
    sal_bridge_lines = {}
    dual_waters = []
    for constituent in line_constituents:
        if (
            constituent[0].split("-")[-1] == "HOH"
            and constituent[2].split("-")[-1] != "HOH"
        ):
            HOH_lines[constituent[0]].append(
                (constituent[2], constituent[3], constituent[4], constituent[5])
            )
        elif (
            constituent[2].split("-")[-1] == "HOH"
            and constituent[0].split("-")[-1] != "HOH"
        ):
            HOH_lines[constituent[2]].append(
                (constituent[0], constituent[1], constituent[4], constituent[5])
            )
        else:
            if constituent[0] not in H_lines:
                H_lines[
                    (
                        constituent[0],
                        constituent[1],
                        constituent[4],
                        constituent[5],
                    )
                ] = [
                    (
                        constituent[2],
                        constituent[3],
                        constituent[4],
                        constituent[5],
                    )
                ]
            else:
                H_lines[
                    (
                        constituent[0],
                        constituent[1],
                        constituent[4],
                        constituent[5],
                    )
                ].append(
                    (
                        constituent[2],
                        constituent[3],
                        constituent[4],
                        constituent[5],
                    )
                )

    fin_HOH_lines = HOH_lines.copy()
    for HOH_line, interactors in HOH_lines.items():
        if len(interactors) < 2:
            del fin_HOH_lines[HOH_line]
    return fin_HOH_lines, H_lines, sal_bridge_lines


def get_parts_from_string(line_part):
    first_0_1 = line_part.find("0")
    res_name1 = line_part[-3:].replace("-", "")  #   B0350-AVAL  C0112ATYR
    # res_num1 = int(line_part[first_0_1:].split("-")[0])
    res_num1 = line_part[1:-3].replace("-", "").lstrip("0")
    chain1 = "".join([i for i in line_part[:-4] if not i.isdigit()])
    return res_name1, res_num1, chain1


def get_individual_constituents(
    line_part1, line_part2, HOH_type: bool = False, H_type: bool = False
):
    final_interactinos = []
    if HOH_type:
        for index, line_part_1 in enumerate(line_part2):
            res_name1, res_num1, chain1 = get_parts_from_string(line_part_1[0])
            for line_part_2 in line_part2[index + 1 :]:
                interaction = {
                    0: {
                        "res_name": None,
                        "res_num": None,
                        "chain": None,
                        "atom": None,
                        "angle": None,
                        "distance": None,
                    },
                    1: {
                        "res_name": None,
                        "res_num": None,
                        "chain": None,
                        "atom": None,
                        "angle": None,
                        "distance": None,
                    },
                }
                res_name2, res_num2, chain2 = get_parts_from_string(line_part_2[0])
                interaction[0]["res_name"] = res_name1
                interaction[0]["res_num"] = res_num1
                interaction[0]["chain"] = chain1
                interaction[0]["atom"] = line_part_1[1]
                interaction[0]["distance"] = line_part_1[2]
                interaction[0]["angle"] = line_part_1[3]
                interaction[1]["res_name"] = res_name2
                interaction[1]["res_num"] = res_num2
                interaction[1]["chain"] = chain2
                interaction[1]["atom"] = line_part_2[1]
                interaction[1]["distance"] = line_part_2[2]
                interaction[1]["angle"] = line_part_2[3]
                final_interactinos.append(interaction)
        return final_interactinos
    elif H_type:
        res_name1, res_num1, chain1 = get_parts_from_string(line_part1[0])

        for part in line_part2:
            interaction = {
                0: {"res_name": None, "res_num": None, "chain": None},
                1: {"res_name": None, "res_num": None, "chain": None},
            }
            res_name2, res_num2, chain2 = get_parts_from_string(part[0])
            interaction[0]["res_name"] = res_name1
            interaction[0]["res_num"] = res_num1
            interaction[0]["chain"] = chain1
            interaction[0]["atom"] = line_part1[1]
            interaction[0]["distance"] = line_part1[2]
            interaction[0]["angle"] = line_part1[3]
            interaction[1]["res_name"] = res_name2
            interaction[1]["res_num"] = res_num2
            interaction[1]["chain"] = chain2
            interaction[1]["atom"] = part[1]
            interaction[1]["distance"] = part[2]
            interaction[1]["angle"] = part[3]
            final_interactinos.append(interaction)
        return final_interactinos


def add_interaction_to_df(
    HOH_lines,
    H_lines,
    chain1,
    chain2,
):
    hydrogen_dict = copy.deepcopy(HYDROGEN_BOND_DICT_STRUCTURE)
    water_dict = copy.deepcopy(WATER_MEDIATED_DICT_STRUCTURE)
    for HOH_line, interactors in HOH_lines.items():
        final_interactions = get_individual_constituents(
            HOH_line, interactors, HOH_type=True
        )
        for final_interaction in final_interactions:
            res1 = final_interaction[0]["res_name"].strip()
            res2 = final_interaction[1]["res_name"].strip()
            res_num1 = final_interaction[0]["res_num"]
            res_num2 = final_interaction[1]["res_num"]
            res1_chain = final_interaction[0]["chain"]
            res2_chain = final_interaction[1]["chain"]
            dist1 = float(final_interaction[0]["distance"])
            dist2 = float(final_interaction[1]["distance"])
            atom1 = final_interaction[0]["atom"]
            atom2 = final_interaction[1]["atom"]
            if (
                res1_chain == chain1
                and res2_chain == chain2
                and dist1 <= WBRIDGE_DIST
                and dist2 <= WBRIDGE_DIST
            ):
                hoh_chain, hoh_res_num = get_items_from_identity(HOH_line.split("-")[0])
                water_dict[RES_NAME_1].append(res1)
                water_dict[RES_NUM_1].append(res_num1)
                water_dict[CHAIN_1].append(res1_chain)
                water_dict[ATOM_1].append(atom1)
                water_dict[RES_NAME_2].append(res2)
                water_dict[RES_NUM_2].append(res_num2)
                water_dict[CHAIN_2].append(res2_chain)
                water_dict[ATOM_2].append(atom2)
                water_dict[DISTANCE_FROM_RES1_NAME].append(dist1)
                water_dict[DISTANCE_FROM_RES2_NAME].append(dist2)
                water_dict[WATER_IDENTITY_NAME].append(f"{hoh_chain}_{hoh_res_num}_HOH")
            elif (
                res2_chain == chain1
                and res1_chain == chain2
                and float(dist1) <= WBRIDGE_DIST
                and float(dist2) <= WBRIDGE_DIST
            ):
                water_dict[RES_NAME_1].append(res2)
                water_dict[RES_NUM_1].append(res_num2)
                water_dict[CHAIN_1].append(res2_chain)
                water_dict[ATOM_1].append(atom2)
                water_dict[RES_NAME_2].append(res1)
                water_dict[RES_NUM_2].append(res_num1)
                water_dict[CHAIN_2].append(res1_chain)
                water_dict[ATOM_2].append(atom1)
                water_dict[DISTANCE_FROM_RES1_NAME].append(dist2)
                water_dict[DISTANCE_FROM_RES2_NAME].append(dist1)
                water_dict[WATER_IDENTITY_NAME].append(HOH_line)

            # self.check_conditions_and_add_interactions(
            #     final_interaction=final_interaction,
            #     interaction_type=INTERACTION_NAMES["Water Mediated"],
            #     HOH_line=HOH_line,
            # )
    for H_line1, interactors in H_lines.items():
        final_interactions = get_individual_constituents(
            H_line1, interactors, H_type=True
        )
        for final_interaction in final_interactions:
            res1 = final_interaction[0]["res_name"].strip()
            res2 = final_interaction[1]["res_name"].strip()
            res_num1 = final_interaction[0]["res_num"]
            res_num2 = final_interaction[1]["res_num"]
            res1_chain = final_interaction[0]["chain"]
            res2_chain = final_interaction[1]["chain"]
            dist1 = float(final_interaction[0]["distance"])
            dist2 = float(final_interaction[1]["distance"])
            atom1 = final_interaction[0]["atom"]
            atom2 = final_interaction[1]["atom"]
            angle = final_interaction[0]["angle"]
            if (
                res2_chain == chain2
                and res1_chain == chain1
                and float(dist2) <= HBOND_DIST
                # and float(angle) >= HBOND_ANGLE
            ):
                hydrogen_dict[RES_NAME_1].append(res1)
                hydrogen_dict[RES_NUM_1].append(res_num1)
                hydrogen_dict[CHAIN_1].append(res1_chain)
                hydrogen_dict[ATOM_1].append(atom1)
                hydrogen_dict[RES_NAME_2].append(res2)
                hydrogen_dict[RES_NUM_2].append(res_num2)
                hydrogen_dict[CHAIN_2].append(res2_chain)
                hydrogen_dict[ATOM_2].append(atom2)
                hydrogen_dict[DISTANCE_NAME].append(dist2)
                hydrogen_dict[DHA_ANGLE_NAME].append(angle)
            elif (
                res2_chain == chain1
                and res1_chain == chain2
                and float(dist2) <= HBOND_DIST
                and float(angle) >= HBOND_ANGLE
            ):
                hydrogen_dict[RES_NAME_1].append(res2)
                hydrogen_dict[RES_NUM_1].append(res_num2)
                hydrogen_dict[CHAIN_1].append(res2_chain)
                hydrogen_dict[ATOM_1].append(atom2)
                hydrogen_dict[RES_NAME_2].append(res1)
                hydrogen_dict[RES_NUM_2].append(res_num1)
                hydrogen_dict[CHAIN_2].append(res1_chain)
                hydrogen_dict[ATOM_2].append(atom1)
                hydrogen_dict[DISTANCE_NAME].append(dist2)
                hydrogen_dict[DHA_ANGLE_NAME].append(angle)
    return hydrogen_dict, water_dict


def get_hydrogen_water_dicts(pdb_file_path, chain1, chain2):
    hb_plus_file = run_hbplus(pdb_file_path)
    hbplus_file_lines = open(hb_plus_file, "r").readlines()
    first_index = None
    all_line_constituents = []
    for index, hbplus_line in enumerate(hbplus_file_lines):
        if hbplus_line == HB_PLUS_FISRT_LINE:
            first_index = index
            break
    for hbplus_line in hbplus_file_lines[first_index + 1 :]:
        line_constituents = [
            hbplus_line[:9].strip(" "),
            hbplus_line[10:13].strip(" "),
            hbplus_line[14:23].strip(" "),
            hbplus_line[24:27].strip(" "),
            hbplus_line[28:32].strip(" "),
            hbplus_line[46:51].strip(" "),
        ]
        all_line_constituents.append(line_constituents)
    HOH_lines, H_lines, salt_bridge_lines = get_hoh_and_h_lines(all_line_constituents)
    hydrogen_dict, water_dict = add_interaction_to_df(
        HOH_lines, H_lines, chain1, chain2
    )
    return hydrogen_dict, water_dict
