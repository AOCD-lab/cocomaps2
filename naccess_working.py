"""
Automates NACCESS-based solvent accessibility calculations for protein–RNA or
protein–protein complexes. Runs NACCESS on full complexes and individual chains,
computes total, polar, and non-polar buried surface areas, and generates detailed
CSV reports for overall and per-chain accessible surface area (ASA) statistics.
"""

import os
import copy

import pandas as pd
from Bio.PDB import NACCESS, PDBParser

from constants import (
    NACCESS_PATH,
    RSA_TABLE,
    NACCESS_FILE_NAME,
    ASA_TABLE_PER_CHAIN,
    RES_NAME_1,
    RES_NUM_1,
    CHAIN_1,
    COMPLEX_ASA_NAME,
    FRESS_ASA_NAME,
    BURRIED_ASA_NAME,
    BURIED_ASA_PERCENTAGE,
    CHAIN1_IDENTIFIER,
    CHAIN2_IDENTIFIER,
)


class Naccess_working:
    def __init__(
        self,
        compound_pdb_file,
        pdb_file_chain1,
        pdb_file_chain2,
        chain1_ress,
        chain2_ress,
        chain_mapping,
    ) -> None:
        self.compound_pdb_path = compound_pdb_file
        self.pdb_chain1_path = pdb_file_chain1
        self.pdb_chain2_path = pdb_file_chain2
        self.chain1_ress = chain1_ress
        self.chain2_ress = chain2_ress
        self.chain_mapping = chain_mapping
        pass

    def run_naccess_on_all(self):
        model_compound = os.path.basename(self.compound_pdb_path)
        model_chain1 = os.path.basename(self.pdb_chain1_path)
        model_chain2 = os.path.basename(self.pdb_chain2_path)
        self.rsa_data_compound, self.asa_data_compound = NACCESS.run_naccess(
            model=model_compound,
            pdb_file=self.compound_pdb_path,
            naccess=NACCESS_PATH,
        )
        self.rsa_data_chain1, self.asa_data_chain1 = NACCESS.run_naccess(
            model=model_chain1,
            pdb_file=self.pdb_chain1_path,
            naccess=NACCESS_PATH,
        )
        self.rsa_data_chain2, self.asa_data_chain2 = NACCESS.run_naccess(
            model=model_chain2,
            pdb_file=self.pdb_chain2_path,
            naccess=NACCESS_PATH,
        )

    def get_area_numbers(self, data):
        text_needed = data[-1].split()
        total, non_polar, polar = (
            float(text_needed[1]),
            float(text_needed[4]),
            float(text_needed[5]),
        )
        return total, non_polar, polar

    def process_area_stats(self):
        self.rsa_dict_compound = NACCESS.process_rsa_data(self.rsa_data_compound)
        self.rsa_dict_chain1 = NACCESS.process_rsa_data(self.rsa_data_chain1)
        self.rsa_dict_chain2 = NACCESS.process_rsa_data(self.rsa_data_chain2)
        (
            self.compound_total_area,
            self.compound_non_polar_area,
            self.compound_polar_area,
        ) = self.get_area_numbers(self.rsa_data_compound)
        self.chain1_total_area, self.chain1_non_polar_area, self.chain1_polar_area = (
            self.get_area_numbers(self.rsa_data_chain1)
        )
        self.chain2_total_area, self.chain2_non_polar_area, self.chain2_polar_area = (
            self.get_area_numbers(self.rsa_data_chain2)
        )
        self.total_burried_area = (
            self.chain1_total_area + self.chain2_total_area - self.compound_total_area
        )
        self.polar_burried_area = (
            self.chain1_polar_area + self.chain2_polar_area - self.compound_polar_area
        )
        self.non_polar_burried_area = (
            self.chain1_non_polar_area
            + self.chain2_non_polar_area
            - self.compound_non_polar_area
        )
        if self.chain1_total_area + self.chain2_total_area == 0:
            self.total_burried_percentage = 100
        else:
            self.total_burried_percentage = (
                self.total_burried_area
                / (self.chain1_total_area + self.chain2_total_area)
                * 100
            )
        if self.total_burried_area == 0:
            self.polar_burried_percentage, self.non_polar_burried_percentage = 100,100
        else:
            self.polar_burried_percentage = (
                self.polar_burried_area
                / (self.total_burried_area)
                * 100
            )
            self.non_polar_burried_percentage = (
                self.non_polar_burried_area
                / (self.total_burried_area)
                * 100
            )

    def make_chain_table(self, chain_num):
        chain_table = copy.deepcopy(ASA_TABLE_PER_CHAIN)
        if chain_num == 1:
            chain = self.chain_mapping[CHAIN1_IDENTIFIER]
            for i in self.rsa_dict_chain1:
                res_num = i[1][1]
                res_name = self.rsa_dict_chain1[i]["res_name"]
                if str(res_num) in self.chain1_ress:
                    free_asa = self.rsa_dict_chain1[i]["all_atoms_abs"]
                    complex_asa = self.rsa_dict_compound[i]["all_atoms_abs"]
                    buried_asa = free_asa - complex_asa
                    if buried_asa >=1:
                        if buried_asa == 0:
                            buried_asa_percentage = 0
                        else:
                            buried_asa_percentage = round((buried_asa / free_asa * 100), 1)
                        chain_table[RES_NAME_1].append(res_name)
                        chain_table[RES_NUM_1].append(res_num)
                        chain_table[CHAIN_1].append(chain)
                        chain_table[COMPLEX_ASA_NAME].append(complex_asa)
                        chain_table[FRESS_ASA_NAME].append(free_asa)
                        chain_table[BURRIED_ASA_NAME].append(round(buried_asa, 2))
                        chain_table[BURIED_ASA_PERCENTAGE].append(buried_asa_percentage)
        if chain_num == 2:
            chain = self.chain_mapping[CHAIN2_IDENTIFIER]
            for i in self.rsa_dict_chain2:
                res_num = i[1][1]
                res_name = self.rsa_dict_chain2[i]["res_name"]
                if str(res_num) in self.chain2_ress:
                    free_asa = self.rsa_dict_chain2[i]["all_atoms_abs"]
                    complex_asa = self.rsa_dict_compound[i]["all_atoms_abs"]
                    buried_asa = free_asa - complex_asa
                    if buried_asa >= 1:
                        if buried_asa == 0:
                            buried_asa_percentage = 0
                        else:
                            buried_asa_percentage = round((buried_asa / free_asa * 100), 1)
                        chain_table[RES_NAME_1].append(res_name)
                        chain_table[RES_NUM_1].append(res_num)
                        chain_table[CHAIN_1].append(chain)
                        chain_table[COMPLEX_ASA_NAME].append(complex_asa)
                        chain_table[FRESS_ASA_NAME].append(free_asa)
                        chain_table[BURRIED_ASA_NAME].append(round(buried_asa, 2))
                        chain_table[BURIED_ASA_PERCENTAGE].append(buried_asa_percentage)
        return chain_table

    def make_table(self):
        rsa_table = copy.deepcopy(RSA_TABLE)
        rsa_table["Value"] = [
            f"{round(self.total_burried_area,2)} / {round((self.total_burried_area/2),2)}",
            round(self.total_burried_percentage, 2),
            f"{round(self.polar_burried_area,2)} / {round((self.polar_burried_area/2),2)}",
            round(self.polar_burried_percentage, 2),
            f"{round(self.non_polar_burried_area,2)} / {round((self.non_polar_burried_area/2),2)}",
            round(self.non_polar_burried_percentage, 2),
        ]

        df = pd.DataFrame(data=rsa_table)
        csv_file_path = f"{self.compound_pdb_path}_{NACCESS_FILE_NAME}"
        df.to_csv(csv_file_path)

    def calculate_stats_and_make_table(self):
        self.run_naccess_on_all()
        self.process_area_stats()
        self.make_table()
        chain1_table = self.make_chain_table(1)
        chain2_table = self.make_chain_table(2)
        csv_file_path1 = f"{self.compound_pdb_path}_ASA_table_chain1.csv"
        csv_file_path2 = f"{self.compound_pdb_path}_ASA_table_chain2.csv"
        df_chain1 = pd.DataFrame(data=chain1_table)
        df_chain2 = pd.DataFrame(data=chain2_table)
        df_chain1.to_csv(csv_file_path1)
        df_chain2.to_csv(csv_file_path2)
