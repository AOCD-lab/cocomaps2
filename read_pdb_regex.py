"""
Parses PDB and mmCIF structure files to extract atomic coordinates and metadata
using regular expressions. Provides utilities for computing interatomic distances
and bond angles, supporting both standard and irregular formatting cases in
structure files. Used as a core input parser for downstream interaction analysis.
"""

import re
import math
import numpy as np
from typing import Optional

from constants import REGEX_NORMAL, REGEX_STRING_NAMES, REGEX_STRING_CHAINS


def disect_line(line):
    raw = line.strip("\n")
    ri = raw[:6].strip()
    raw_n = raw[6:]
    matched = re.search(
        r"(\d+)\s*([A-Za-z0-9]+\'?)\s+([A-Za-z0-9]{1,4})\s*([A-Z0-9])\s*(\w*\d+\w*)\s+(-?\d+.?\d+)\s*(-?\d+.?\d+)\s*(-?\d+.?\d+)\s+(\d+.?\d+)\s*(\d+.?\d+)\s+(\w+)",
        raw_n,
    )
    if matched and ri in ["ATOM", "HETATM"]:
        p = [
            ri,
            matched.group(1),
            matched.group(2),
            matched.group(3),
            matched.group(4),
            matched.group(5),
            float(matched.group(6)),
            float(matched.group(7)),
            float(matched.group(8)),
            float(matched.group(9)),
            matched.group(10),
            matched.group(11),
        ]
        return p


def read_pdb(
    file: str,
    cif_bool: Optional[bool] = False,
    all_res_names: Optional[list] = None,
    all_chain_names: Optional[list] = None,
):
    lines = []
    parsed = []
    if not cif_bool:
        with open(file) as fh:
            for line in fh.readlines():
                raw = line.strip("\n")
                ri = raw[:6].strip()
                raw_n = raw[6:]
                matched = re.search(
                    # r"(\d+)\s*([A-Za-z0-9]+\'?)\s+([A-Za-z0-9]{1,4})\s*([a-zA-Z0-9])\s*(\w*\d+\w*)\s+(-?\d+.?\d+)\s*(-?\d+.?\d+)\s*(-?\d+.?\d+)\s+(\d+.?\d+)\s*(\d+.?\d+)\s+(\w+)",
                    r"(\d+)\s*([A-Za-z0-9]+\'{0,2}\'?)\s+([A-Za-z0-9]{1,4})\s*([a-zA-Z0-9@])\s*(\w*\d+\w*)\s+(-?\d+\.?\d+)\s*(-?\d+\.?\d+)\s*(-?\d+\.?\d+)\s+(\d+\.?\d+)\s*(\d+\.?\d+)\s+(\w+)",
                    raw_n,
                )
                if matched and ri in ["ATOM", "HETATM"]:
                    try:
                        p = [
                            ri,
                            matched.group(1),
                            matched.group(2),
                            matched.group(3),
                            matched.group(4),
                            matched.group(5),
                            float(line[30:38].strip()),  # Extract X coordinate
                            float(line[38:46].strip()),  # Extract Y coordinate
                            float(line[46:54].strip()),
                            float(matched.group(9)),
                            matched.group(10),
                            matched.group(11),
                        ]
                        # print(p)
                        parsed.append(p)
                        lines.append(raw)
                    except Exception as e:
                        print(f"reading the pdb error {e}")
    else:
        with open(file) as fh:
            for line in fh.readlines():
                raw = line.strip("\n")
                ri = raw[:6].strip()
                raw_n = raw[6:]
                matched = re.search(
                    # r"(\d+)\s*([A-Za-z0-9]+\'?)\s+([A-Za-z0-9]{1,4})\s*([a-zA-Z0-9])\s*(\w*\d+\w*)\s+(-?\d+.?\d+)\s*(-?\d+.?\d+)\s*(-?\d+.?\d+)\s+(\d+.?\d+)\s*(\d+.?\d+)\s+(\w+)",
                    REGEX_NORMAL,
                    raw_n,
                )
                res_name_iterator = 4
                while (
                    matched
                    and matched.group(3) not in all_res_names
                    and res_name_iterator > 0
                ):
                    try:
                        new_regex_string = REGEX_STRING_NAMES % (res_name_iterator)
                        res_name_iterator -= 1
                        matched = re.search(
                            rf"{new_regex_string}",
                            raw_n,
                        )
                    except:
                        matched = None
                chain_iterator = 1
                while (
                    matched
                    and chain_iterator < 3
                    and matched.group(4) not in all_chain_names
                ):
                    try:
                        new_regex_string = REGEX_STRING_CHAINS % (
                            res_name_iterator,
                            chain_iterator,
                        )
                        chain_iterator += 1
                        matched = re.search(
                            rf"{new_regex_string}",
                            raw_n,
                        )
                    except:
                        matched = None

                if matched and ri in ["ATOM", "HETATM"]:
                    try:
                        p = [
                            ri,
                            matched.group(1),
                            matched.group(2),
                            matched.group(3),
                            matched.group(4),
                            matched.group(5),
                            float(line[30:39].strip()),  # Extract X coordinate
                            float(line[39:46].strip()),  # Extract Y coordinate
                            float(line[46:54].strip()),
                            float(matched.group(9)),
                            matched.group(10),
                            matched.group(11),
                        ]
                        parsed.append(p)
                        lines.append(raw)
                    except Exception as e:
                        print(f"reading the pdb error {e}")

    return parsed, lines


def distance(p1, p2):
    x = p1[6] - p2[6]
    y = p1[7] - p2[7]
    z = p1[8] - p2[8]
    d = math.sqrt(x * x + y * y + z * z)
    return d


def get_angle(line1, line2, line3):
    a = np.array([line1[6], line1[7], line1[8]])
    b = np.array([line2[6], line2[7], line2[8]])
    c = np.array([line3[6], line3[7], line3[8]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


# def main():
#     for file in os.listdir('/Users/kanavkalra/Desktop/all_pdb/'):
#      #for file in ['3u5h.pdb']:
#         if file.endswith('.pdb'):
#             print(file)
#             parsed, lines = read_pdb('/home/kalrak/all_pdb/' + file)
#             print(parsed)


# if __name__ == '__main__':
#     main()

# read_pdb(
#     "/Users/utkarshkalra/kaust/clear/new_cmd_interactions/300pdbs/RQH-flip-hydrogen/6ey6.pdb-flip.pdb"
# )
