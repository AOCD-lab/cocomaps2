"""
Detects disulfide (S–S) and salt bridge interactions within molecular structures
based on residue and atom types defined in residue-specific rules. Applies
distance cutoffs to identify valid covalent and electrostatic contacts, and
stores the results in structured dictionaries for further analysis.
"""

import copy

from residues import S_BOND, SALT_BRIDGE_CONDITION1, SALT_BRIDGE_CONDITION2
from constants import (
    SSBOND_DIST,
    SS_BOND_DICT_STRUCTURE,
    SBRIDGE_DIST,
    SALT_BRIDGE_DICT_STRUCTURE,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_NAME,
)

from update_constants import tc

SBRIDGE_DIST = tc.SBRIDGE_DIST
SSBOND_DIST = tc.SSBOND_DIST


def check_ss_bond(res_1, atom1, res_2, atom2, distance):
    if (
        res_1 in S_BOND.keys()
        and res_2 in S_BOND.keys()
        and atom1 in S_BOND.values()
        and atom2 in S_BOND.values()
        and distance <= SSBOND_DIST
    ):
        return True
    return False


def check_salt_bridge(res_1, atom1, res_2, atom2, distance):
    if distance <= SBRIDGE_DIST:
        if (
            res_1 in SALT_BRIDGE_CONDITION1["side1"]
            and res_2 in SALT_BRIDGE_CONDITION1["side2"]
            and atom1 in SALT_BRIDGE_CONDITION1["side1"][res_1]
            and atom2 in SALT_BRIDGE_CONDITION1["side2"][res_2]
        ):
            return True
        elif (
            res_2 in SALT_BRIDGE_CONDITION1["side1"]
            and res_1 in SALT_BRIDGE_CONDITION1["side2"]
            and atom2 in SALT_BRIDGE_CONDITION1["side1"][res_2]
            and atom1 in SALT_BRIDGE_CONDITION1["side2"][res_1]
        ):
            return True
        elif (
            res_1 in SALT_BRIDGE_CONDITION2["side1"]
            and atom1 in SALT_BRIDGE_CONDITION2["side1"][res_1]
            and res_2 in SALT_BRIDGE_CONDITION2["side2"]
            and atom2 in SALT_BRIDGE_CONDITION2["allowed_atoms_side2"]
        ):
            return True
        elif (
            res_2 in SALT_BRIDGE_CONDITION2["side1"]
            and atom2 in SALT_BRIDGE_CONDITION2["side1"][res_2]
            and res_1 in SALT_BRIDGE_CONDITION2["side2"]
            and atom1 in SALT_BRIDGE_CONDITION2["allowed_atoms_side2"]
        ):
            return True
    return False


def get_ss_bond_and_sbridge_interactions(full_df):
    ss_bond_interactions = copy.deepcopy(SS_BOND_DICT_STRUCTURE)
    salt_bridge_interactions = copy.deepcopy(SALT_BRIDGE_DICT_STRUCTURE)
    for index, row in full_df.iterrows():
        res_1 = row[RES_NAME_1]
        atom1 = row[ATOM_1]
        atom2 = row[ATOM_2]
        res_num_1 = row[RES_NUM_1]
        chain_1 = row[CHAIN_1]
        res_2 = row[RES_NAME_2]
        res_num_2 = row[RES_NUM_2]
        chain_2 = row[CHAIN_2]
        distance = row[DISTANCE_NAME]
        if check_ss_bond(res_1, atom1, res_2, atom2, distance):
            ss_bond_interactions[RES_NAME_1].append(res_1)
            ss_bond_interactions[RES_NUM_1].append(res_num_1)
            ss_bond_interactions[CHAIN_1].append(chain_1)
            ss_bond_interactions[ATOM_1].append(atom1)
            ss_bond_interactions[RES_NAME_2].append(res_2)
            ss_bond_interactions[RES_NUM_2].append(res_num_2)
            ss_bond_interactions[CHAIN_2].append(chain_2)
            ss_bond_interactions[ATOM_2].append(atom2)
            ss_bond_interactions[DISTANCE_NAME] = distance
        if check_salt_bridge(res_1, atom1, res_2, atom2, distance):
            salt_bridge_interactions[RES_NAME_1].append(res_1)
            salt_bridge_interactions[RES_NUM_1].append(res_num_1)
            salt_bridge_interactions[CHAIN_1].append(chain_1)
            salt_bridge_interactions[ATOM_1].append(atom1)
            salt_bridge_interactions[RES_NAME_2].append(res_2)
            salt_bridge_interactions[RES_NUM_2].append(res_num_2)
            salt_bridge_interactions[CHAIN_2].append(chain_2)
            salt_bridge_interactions[ATOM_2].append(atom2)
            salt_bridge_interactions[DISTANCE_NAME].append(distance)

    return ss_bond_interactions, salt_bridge_interactions
