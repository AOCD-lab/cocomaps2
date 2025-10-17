"""
Main processing pipeline for analyzing biomolecular interactions between two chains.

This module orchestrates the full workflow — from reading and preprocessing PDB data
to detecting, classifying, and recording all major non-covalent and covalent interactions.
It integrates all specialized submodules (π–π, hydrogen bonds, salt bridges, van der Waals,
metal-mediated, and accessibility analyses) into a unified data processing sequence.

Key operations include:
- Constructing atom-pair distance tables and identifying ring systems.
- Detecting specific interaction types (H-bonds, halogen bonds, π-interactions, etc.)
  based on geometric and chemical criteria.
- Running NACCESS to compute buried surface areas and solvent accessibility metrics.
- Generating individual CSV reports per interaction type and a consolidated summary file.
- Producing formatted, truncated PDB structures containing only relevant atoms.

The `process()` function serves as the main entry point and coordinates the entire
interaction detection and output generation pipeline for a given protein/RNA complex.
"""

import os
from typing import Optional

from utils import (
    initial_filter_data_frame_creator,
    make_truncated_pdb,
    remove_salt_duplicates,
    filter_hydrogen_salt_dicts,
)

from pi_fucntions import Ring_params
from pi_interactions import get_all_pi_interactions, get_pi_pi_interactions
from metal_mediated import get_intercations_metal
from c_h_pi_interactions import get_c_h_pi_interactions
from n_h_pi_interactions import get_nso_h_pi_interactions
from weak_h_bond_interactions import get_weak_h_bond_interactions
from halogen_interactions import get_halogen_interactions
from hb_plus import get_hydrogen_water_dicts
from polar_apolar_vdw_interactions import get_apolar_vdw, get_polar_vdw
from ss_bond_sbridge_interactions import get_ss_bond_and_sbridge_interactions
from df_works import DF_Works
from naccess_working import Naccess_working
from ligand_identifier import identify_ligands

from constants import DISTANCE_NAME, ATOM_1, ATOM_2, ATOM_1_ORG_NAME, ATOM_2_ORG_NAME

from utils import get_c_h_n_h_lines

from constants import (
    CHAIN1_IDENTIFIER,
    CHAIN2_IDENTIFIER,
    WEAK_H_BONDS_NAME,
    WEAK_H_BONDS_PRINT_NAME,
    METAL_MEDIATED_NAME,
    METAL_MEDIATED_PRINT_NAME,
    PI_PI_NAME,
    PI_PI_PRINT_NAME,
    LONE_PAIR_PI_NAME,
    LONE_PAIR_PI_PRINT_NAME,
    ANION_PI_NAME,
    ANION_PI_PRINT_NAME,
    CATION_PI_NAME,
    CATION_PI_PRINT_NAME,
    AMINO_PI_NAME,
    AMINO_PI_PRINT_NAME,
    C_H_PI_NAME,
    C_H_PI_PRINT_NAME,
    N_H_PI_NAME,
    N_H_PI_PRINT_NAME,
    APOLAR_NAME,
    APOLAR_PRINT_NAME,
    POLAR_NAME,
    POLAR_PRINT_NAME,
    HYDROGEN_NAME,
    HYDROGEN_PRINT_NAME,
    HALOGEN_NAME,
    HALOGEN_PRINT_NAME,
    WATER_MEDIATED_NAME,
    WATER_MEDIATED_PRINT_NAME,
    SALT_BRIDGE_NAME,
    SALT_BRIDGE_PRINT_NAME,
    SS_BOND_NAME,
    SS_BOND_PRINT_NAME,
    RES_NUM_1,
    RES_NUM_2,
    WATER_IDENTITY_NAME,
    METAL_IDENTITY_NAME,
)


def process(
    org_pdb_path,
    org_chains,
    chains_lines: dict,
    chains_mapping: dict,
    all_pdb_parsed_lines: list,
    all_pdb_formatted_file_path: str,
    compound_pdb_path: str,
    chain1_pdb_path: str,
    chain2_pdb_path: str,
):

    ligands_identified = identify_ligands(org_pdb_path)
    chains_parsed_lines = []
    for chain, chain_lines in chains_lines.items():
        chains_parsed_lines += chain_lines

    full_df, df, rings, metal_oriented_lines = initial_filter_data_frame_creator(
        chains_lines
    )
    ring_params = Ring_params()
    ring_params.make_chains_ring_comps(rings)
    pi_pi_interactions = get_pi_pi_interactions(ring_params, df)
    (
        lone_pair_pi_interactions,
        anion_pi_interactions,
        cation_pi_interactions,
        amino_pi_interactions,
    ) = get_all_pi_interactions(full_df, ring_params, df)
    metal_mediated_interactions = get_intercations_metal(
        all_pdb_parsed_lines, chains_mapping, metal_oriented_lines
    )

    c_h_combos, n_h_combos = get_c_h_n_h_lines(chains_lines)
    c_h_pi_interactions = get_c_h_pi_interactions(df, c_h_combos, ring_params)
    nso_h_pi_interactions = get_nso_h_pi_interactions(df, n_h_combos, ring_params)
    weak_h_bond_interactions = get_weak_h_bond_interactions(
        chains_parsed_lines, [CHAIN1_IDENTIFIER, CHAIN2_IDENTIFIER]
    )
    halogen_interactions = get_halogen_interactions(
        chains_parsed_lines, [CHAIN1_IDENTIFIER, CHAIN2_IDENTIFIER]
    )
    hydrogen_bond_interactions, water_mediated_interactions = get_hydrogen_water_dicts(
        all_pdb_formatted_file_path, CHAIN1_IDENTIFIER, CHAIN2_IDENTIFIER
    )
    apolar_vdw_interactions = get_apolar_vdw(full_df, pi_pi_interactions)
    polar_vdw_interactions = get_polar_vdw(
        full_df,
        pi_pi_interactions,
        hydrogen_bond_interactions,
        halogen_interactions,
        weak_h_bond_interactions,
    )
    ss_bond_interactions, salt_bridge_interactions = (
        get_ss_bond_and_sbridge_interactions(full_df)
    )
    hydrogen_bond_interactions = filter_hydrogen_salt_dicts(
        hydrogen_bond_interactions, salt_bridge_interactions
    )
    salt_bridge_interactions = remove_salt_duplicates(salt_bridge_interactions)
    try:
        naccess_working = Naccess_working(
            compound_pdb_file=compound_pdb_path,
            pdb_file_chain1=chain1_pdb_path,
            pdb_file_chain2=chain2_pdb_path,
            chain1_ress=df[RES_NUM_1].tolist(),
            chain2_ress=df[RES_NUM_2].tolist(),
            chain_mapping=chains_mapping,
        )
        naccess_working.calculate_stats_and_make_table()
    except Exception as e:
        print(f"[WARNING] Naccess error: {e}")

    df_works = DF_Works(org_pdb_path, df, org_chains)
    df_works.final_df_reverse_chain_mapping(chains_mapping)
    sb_df = df_works.save_inidividual_csv(
        ss_bond_interactions,
        SS_BOND_NAME,
        0,
        chains_mapping,
    )
    df_works.add_to_final_df(sb_df, SS_BOND_PRINT_NAME)

    s_df = df_works.save_inidividual_csv(
        salt_bridge_interactions, SALT_BRIDGE_NAME, 1, chains_mapping
    )
    df_works.add_to_final_df(s_df, SALT_BRIDGE_PRINT_NAME)

    hb_df = df_works.save_inidividual_csv(
        hydrogen_bond_interactions, HYDROGEN_NAME, 2, chains_mapping
    )
    df_works.add_to_final_df(hb_df, HYDROGEN_PRINT_NAME)

    w_df = df_works.save_inidividual_csv(
        water_mediated_interactions, WATER_MEDIATED_NAME, 3, chains_mapping
    )
    df_works.add_to_final_df(w_df, WATER_MEDIATED_PRINT_NAME)

    h_df = df_works.save_inidividual_csv(
        weak_h_bond_interactions, WEAK_H_BONDS_NAME, 4, chains_mapping
    )
    df_works.add_to_final_df(h_df, WEAK_H_BONDS_PRINT_NAME)

    halo_df = df_works.save_inidividual_csv(
        halogen_interactions, HALOGEN_NAME, 5, chains_mapping
    )
    df_works.add_to_final_df(halo_df, HALOGEN_PRINT_NAME)

    m_df = df_works.save_inidividual_csv(
        metal_mediated_interactions, METAL_MEDIATED_NAME, 6, chains_mapping
    )
    df_works.add_to_final_df(m_df, METAL_MEDIATED_PRINT_NAME)

    p_df = df_works.save_inidividual_csv(
        pi_pi_interactions, PI_PI_NAME, 7, chains_mapping
    )
    df_works.add_to_final_df(p_df, PI_PI_PRINT_NAME)

    lp_df = df_works.save_inidividual_csv(
        lone_pair_pi_interactions, LONE_PAIR_PI_NAME, 8, chains_mapping
    )
    df_works.add_to_final_df(lp_df, LONE_PAIR_PI_PRINT_NAME)

    ap_df = df_works.save_inidividual_csv(
        anion_pi_interactions, ANION_PI_NAME, 9, chains_mapping
    )
    df_works.add_to_final_df(ap_df, ANION_PI_PRINT_NAME)

    cp_df = df_works.save_inidividual_csv(
        cation_pi_interactions, CATION_PI_NAME, 10, chains_mapping
    )
    df_works.add_to_final_df(cp_df, CATION_PI_PRINT_NAME)

    amino_pi_df = df_works.save_inidividual_csv(
        amino_pi_interactions, AMINO_PI_NAME, 11, chains_mapping
    )
    df_works.add_to_final_df(amino_pi_df, AMINO_PI_PRINT_NAME)

    nhp_df = df_works.save_inidividual_csv(
        nso_h_pi_interactions, N_H_PI_NAME, 12, chains_mapping
    )
    df_works.add_to_final_df(nhp_df, N_H_PI_PRINT_NAME)

    chp_df = df_works.save_inidividual_csv(
        c_h_pi_interactions, C_H_PI_NAME, 13, chains_mapping
    )
    df_works.add_to_final_df(chp_df, C_H_PI_PRINT_NAME)

    pl_df = df_works.save_inidividual_csv(
        polar_vdw_interactions, POLAR_NAME, 14, chains_mapping
    )
    df_works.add_to_final_df(pl_df, POLAR_PRINT_NAME)

    al_df = df_works.save_inidividual_csv(
        apolar_vdw_interactions, APOLAR_NAME, 15, chains_mapping
    )
    df_works.add_to_final_df(al_df, APOLAR_PRINT_NAME)

    df_works.add_clash_to_final_df()
    df_works.add_proximal_to_final_df()
    df_works.save_clash_proximal_files()
    df_works.drop_unwanted_rows(ligands_identified)
    df_works.drop_columns_from_final_df(
        [DISTANCE_NAME, ATOM_1, ATOM_2, ATOM_1_ORG_NAME, ATOM_2_ORG_NAME]
    )
    df_works.make_summary_table()
    df_works.save_final_df()

    make_truncated_pdb(
        water_mediated_interactions[WATER_IDENTITY_NAME],
        metal_mediated_interactions[METAL_IDENTITY_NAME],
        chains_mapping,
        all_pdb_parsed_lines,
        compound_pdb_path,
    )
