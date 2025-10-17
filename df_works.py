"""
This module manages all DataFrame-related operations for molecular interaction data.
It handles combining, labeling, and saving interaction results (e.g., clashes, proximals),
updates summary tables, filters unwanted entries like water and metals, and generates
final CSV reports summarizing detected residue–residue interactions between chains.
"""

import os
import copy
import pandas as pd
from typing import Optional
from bidict import bidict
from constants import (
    TYPE_COLUMN_NAME,
    PROXIMAL_PRINT_NAME,
    BASIC_INTERACTION_DICT_STRUCTURE,
    CLASH_PRINT_NAME,
    PROXIMAL_PRINT_NAME,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    ATOM_1,
    RES_NAME_2,
    RES_NUM_2,
    CHAIN_2,
    ATOM_2,
    DISTANCE_NAME,
    PROXIMAL_NAME,
    CLASH_NAME,
    INTERNAL_CLASH_CHECK_INTERCATION_TYPES,
    SUMMARY_TABLE,
    SMALL_SUMMARY_TABLE,
    DROP_ROWS
)
from residues import RADII, METALS
from utils import reverse_mapping_chains_in_df, check_clash_in_row


class DF_Works:
    def __init__(
        self,
        pdb_path,
        final_df,
        chains_org,
    ) -> None:
        self.pdb_path = pdb_path
        self.base_name = f"{pdb_path}_{chains_org[0]}_{chains_org[1]}"
        self.final_df = final_df
        self.final_df = self.final_df.drop(["x1", "y1", "z1", "x2", "y2", "z2"], axis=1)
        self.final_df[TYPE_COLUMN_NAME] = ""
        self.proximal_df = pd.DataFrame(
            columns=list(BASIC_INTERACTION_DICT_STRUCTURE.keys())
        )
        self.clash_df = pd.DataFrame(
            columns=list(BASIC_INTERACTION_DICT_STRUCTURE.keys())
        )
        self.summary_table = copy.deepcopy(SUMMARY_TABLE)
        self.small_summary_table = copy.deepcopy(SMALL_SUMMARY_TABLE)
        pass

    def save_inidividual_csv(
        self,
        interaction_dict,
        interaction_type,
        ind,
        chain_mapping: Optional[bidict] = None,
    ):
        csv_name = f"{self.base_name}_{interaction_type}.csv"
        temp_df = pd.DataFrame(data=interaction_dict)
        temp_df = temp_df.drop_duplicates()
        temp_df = reverse_mapping_chains_in_df(chain_mapping, temp_df, interaction_type)
        indices_to_delete = []
        if interaction_type in INTERNAL_CLASH_CHECK_INTERCATION_TYPES:
            for index, row in temp_df.iterrows():
                if ((self.final_df[RES_NUM_1] == row[RES_NUM_1]) & (self.final_df[RES_NUM_2] == row[RES_NUM_2])).any():
                    if check_clash_in_row(row):
                        temp_df.loc[index, DISTANCE_NAME] = f"{row[DISTANCE_NAME]} *"
                else:
                    indices_to_delete.append(index)
        temp_df.drop(index=indices_to_delete, inplace=True)
        temp_df = temp_df[~temp_df[RES_NAME_1].isin(DROP_ROWS)]
        temp_df = temp_df[~temp_df[RES_NAME_2].isin(DROP_ROWS)]
        temp_df.to_csv(csv_name)
        self.summary_table["Value"][ind] = int(temp_df.shape[0])
        return temp_df

    def add_to_final_df(self, interaction_df, interaction_type):
        # inetractions_len = len(interaction_dict["RNA Res. Name"])
        for index, row in interaction_df.iterrows():
            res_name_1 = row[RES_NAME_1]
            res_num_1 = row[RES_NUM_1]
            chain_1 = row[CHAIN_1]
            res_name_2 = row[RES_NAME_2]
            res_num_2 = row[RES_NUM_2]
            chain_2 = row[CHAIN_2]
            condition = (
                (self.final_df[RES_NAME_1] == res_name_1)
                & (self.final_df[RES_NUM_1] == res_num_1)
                & (self.final_df[CHAIN_1] == chain_1)
                & (self.final_df[RES_NAME_2] == res_name_2)
                & (self.final_df[RES_NUM_2] == res_num_2)
                & (self.final_df[CHAIN_2] == chain_2)
            )

            # Find the index of the row meeting the condition
            index = self.final_df.loc[condition].index
            if not index.any():
                new_row = {
                    RES_NAME_1: res_name_1,
                    RES_NUM_1: res_num_1,
                    CHAIN_1: chain_1,
                    ATOM_1: "-",
                    RES_NAME_2: res_name_2,
                    RES_NUM_2: res_num_2,
                    CHAIN_2: chain_2,
                    ATOM_2: "-",
                    DISTANCE_NAME: 0,  # because this will be removed anyway so lets not get into more error traps
                    TYPE_COLUMN_NAME: f"; {interaction_type}",
                }
                self.final_df.loc[len(self.final_df)] = new_row
            # Append value to a specific column in that row
            column_to_append = TYPE_COLUMN_NAME
            if (
                interaction_type in INTERNAL_CLASH_CHECK_INTERCATION_TYPES
                and "*" in str(row[DISTANCE_NAME])
            ):
                value_to_append = (
                    self.final_df.loc[index, column_to_append]
                    + f"; {interaction_type}*"
                )  # Concatenate with "black"
            else:
                value_to_append = (
                    self.final_df.loc[index, column_to_append] + f"; {interaction_type}"
                )  # Concatenate with "black"
            self.final_df.loc[index, column_to_append] = value_to_append

    def add_proximal_to_final_df(self):
        for index, row in self.final_df.iterrows():
            if not row[TYPE_COLUMN_NAME]:
                self.final_df.loc[index, TYPE_COLUMN_NAME] = PROXIMAL_PRINT_NAME
                self.proximal_df = pd.concat([self.proximal_df, pd.DataFrame([row])])

    def final_df_reverse_chain_mapping(self, chain_mapping):
        self.final_df = reverse_mapping_chains_in_df(chain_mapping, self.final_df)

    def save_final_df(self):
        file_name = f"{self.base_name}_final_file.csv"
        self.final_df.to_csv(file_name)

    def add_clash_to_final_df(self):
        for index, row in self.final_df.iterrows():
            if not row[TYPE_COLUMN_NAME]:
                if check_clash_in_row(row):
                    self.final_df.loc[index, TYPE_COLUMN_NAME] = CLASH_PRINT_NAME
                    self.clash_df = pd.concat([self.clash_df, pd.DataFrame([row])])

    def drop_columns_from_final_df(self, column_names):
        self.final_df = self.final_df.drop(columns=column_names)

    def make_summary_table(self):
        chain1_res_num = len(set(self.final_df[RES_NUM_1].tolist()))
        chain2_res_num = len(set(self.final_df[RES_NUM_2].tolist()))
        self.small_summary_table["Value"][0] = chain1_res_num + chain2_res_num
        self.small_summary_table["Value"][1] = chain1_res_num
        self.small_summary_table["Value"][2] = chain2_res_num
        small_summary_df = pd.DataFrame(data=self.small_summary_table)
        small_summary_df.to_csv(f"{self.base_name}_small_summary.csv")
        summary_df = pd.DataFrame(data=self.summary_table)
        summary_df.to_csv(f"{self.base_name}_summary_table.csv")

    def drop_unwanted_rows(self, ligands_identified):
        drop_index = []
        for index, row in self.final_df.iterrows():
            if (row[RES_NAME_1], row[CHAIN_1], row[RES_NUM_1]) in ligands_identified:
                drop_index.append(index)
            elif (row[RES_NAME_2], row[CHAIN_2], row[RES_NUM_2]) in ligands_identified:
                drop_index.append(index)
        self.final_df.drop(drop_index, inplace=True)


        self.final_df = self.final_df[~self.final_df[RES_NAME_1].isin(DROP_ROWS)]
        self.final_df = self.final_df[~self.final_df[RES_NAME_2].isin(DROP_ROWS)]
    
        for met in METALS:
            self.final_df = self.final_df[self.final_df[RES_NAME_1] != met]
            self.final_df = self.final_df[self.final_df[RES_NAME_2] != met]
        

    def save_clash_proximal_files(self):
        csv_proximal_name = f"{self.base_name}_{PROXIMAL_NAME}.csv"
        csv_clash_name = f"{self.base_name}_{CLASH_NAME}.csv"
        self.proximal_df.to_csv(csv_proximal_name)
        self.clash_df.to_csv(csv_clash_name)
        self.summary_table["Value"][-2] = len(self.clash_df)
        self.summary_table["Value"][-1] = len(self.proximal_df)
