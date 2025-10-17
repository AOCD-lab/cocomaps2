"""
Defines the Threshhold_constants class, which stores and manages all configurable
interaction cutoff and angle parameters used across the pipeline. It initializes
values from `constants.py`, allows runtime updates via JSON input, and provides
a flexible interface for modifying thresholds without editing source code.
"""

import json
import constants
import inspect


class Threshhold_constants:
    def __init__(self) -> None:
        self.pdb_file = constants.pdb_file
        self.REDUCE_BOOL = constants.REDUCE_BOOL
        self.chains_set_1 = constants.chains_set_1
        self.chains_set_2 = constants.chains_set_2
        self.ranges_1 = constants.ranges_1
        self.ranges_2 = constants.ranges_2
        self.CUT_OFF = constants.CUT_OFF
        self.HBOND_DIST = constants.HBOND_DIST
        self.HBOND_ANGLE = constants.HBOND_ANGLE
        self.SBRIDGE_DIST = constants.SBRIDGE_DIST
        self.WBRIDGE_DIST = constants.WBRIDGE_DIST
        self.SSBOND_DIST = constants.SSBOND_DIST
        self.HALOGEN_THETA1 = constants.HALOGEN_THETA1
        self.HALOGEN_THETA2 = constants.HALOGEN_THETA2
        self.METAL_DIST = constants.METAL_DIST
        self.PI_PI_DIST = constants.PI_PI_DIST
        self.PI_PI_THETA = constants.PI_PI_THETA
        self.PI_PI_GAMMA = constants.PI_PI_GAMMA
        self.ANION_PI_DIST = constants.ANION_PI_DIST
        self.LONEPAIR_PI_DIST = constants.LONEPAIR_PI_DIST
        self.CATION_PI_DIST = constants.CATION_PI_DIST
        self.POLAR_PI_DIST = constants.AMINO_PI_DIST
        self.CH_ON_DIST = constants.CH_ON_DIST
        self.CH_ON_ANGLE = constants.CH_ON_ANGLE
        self.C_H_PI_DIST = constants.C_H_PI_DIST
        self.C_H_PI_THETA1 = constants.C_H_PI_THETA1
        self.C_H_PI_THETA2 = constants.C_H_PI_THETA2
        self.NSOH_PI_THETA1 = constants.NSOH_PI_THETA1
        self.NSOH_PI_THETA2 = constants.NSOH_PI_THETA2
        self.NSOH_PI_DIST = constants.NSOH_PI_DIST
        self.APOLAR_TOLERANCE = constants.APOLAR_TOLERANCE
        self.POLAR_TOLERANCE = constants.POLAR_TOLERANCE
        self.SBRIDGE_DIST = constants.SBRIDGE_DIST
        pass

    def get_var_name(var):
        callers_local_vars = inspect.currentframe().f_back.f_locals
        return [
            var_name
            for var_name, var_val in callers_local_vars.items()
            if var_val is var
        ][0]

    def update_variables(self, input_dict):
        for var_name, new_value in input_dict.items():
            if hasattr(self, var_name):
                setattr(self, var_name, new_value)
                print(f"Variable '{var_name}' updated to {new_value}.")
            else:
                print(f"Variable '{var_name}' not found.")

    def read_json_file(self, json_file: str):
        fl = open(json_file)
        data = json.load(fl)
        self.update_variables(data)


tc = Threshhold_constants()
