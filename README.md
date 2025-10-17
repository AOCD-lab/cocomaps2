# coco_2



## ⚙️ Installation


### 1. Create and activate a Conda environment

```bash
conda create -n coco2 python=3.12
conda activate coco2
```

### 2. Install Python dependencies

```bash
conda config --append channels conda-forge
conda install --file coco2/conda_requirements.txt --yes
```

### 3. External tools

| Tool | Purpose |
|------|----------|
| **HBPLUS** | Hydrogen bond detection | 
| **HBADD** | Hydrogen bond information to be used by HBPLUS |
| **NACCESS** | Accessible surface area | 
| **REDUCE** | Hydrogen optimization and Addition | 


All the links to these softwares can be found in the references.
### 4. Paths set up

Executable paths of the sofwares mentioned above should be entered to specific variables in `constants.py` file.

```python
HB_PLUS_PATH = "path_to_hblpus"
HBADD_PATH = "path_to_hbadd"
HBADD_EMPTY_FILE_PATH = "path_to_an_empty_.txt_file"
REDUCE_COMMAND = "path_to_reduce -FLIP %s > %s"
NACCESS_PATH = "path_to_naccess"
 ``` 



---



## ▶️ Usage

To run COCO2 you need to set up a .json input file:

The File format:
```json
{
    "pdb_file": "path to the input pdb/cif file",
    "chains_set_1": [
        chain_set 1
    ],
    "chains_set_2": [
        chain_set 2
    ],
    "ranges_1": [
        [
            start_res_number for chain_set_1[0],
            end_res_number for chain_set1_[0]
        ]
    ],
    "ranges_2": [
        [
            start_res_number for chain_set_2[0],
            end_res_number for chain_set_2[0]
        ],
        [
            start_res_number for chain_set_2[1],
            end_res_number for chain_set_2[1]
        ]
    ],
    "HBOND_DIST": 3.9,
    "HBOND_ANGLE": 90,
    "SBRIDGE_DIST": 4.5,
    "WBRIDGE_DIST": 3.9,
    "CH_ON_DIST": 3.6,
    "CH_ON_ANGLE": 110,
    "CUT_OFF": 5,
    "APOLAR_TOLERANCE": 0.5,
    "POLAR_TOLERANCE": 0.5,
    "PI_PI_DIST": 5.5,
    "PI_PI_THETA": 80,
    "PI_PI_GAMMA": 90,
    "ANION_PI_DIST": 5,
    "LONEPAIR_PI_DIST": 5,
    "AMINO_PI_DIST": 5,
    "CATION_PI_DIST": 5,
    "METAL_DIST": 3.2,
    "HALOGEN_THETA1": 165,
    "HALOGEN_THETA2": 120,
    "C_H_PI_DIST": 5.0,
    "C_H_PI_THETA1": 120,
    "C_H_PI_THETA2": 30,
    "NSOH_PI_DIST": 4.5,
    "NSOH_PI_THETA1": 120,
    "NSOH_PI_THETA2": 30
}
```
parameter values can be changed as desired


To run
```bash
python begin.py input_json_path.json
```


## 📊 Output Files

After running, COCO2 generates multiple `.csv` files in the output directory:
'*' represents the specific prefix name to the files according to the input pdb name and chains
| File Name | Description |
|------------|-------------|
| `*_final_file.csv` | All contacts detected |
| `*_Amino_pi.csv` | All Amino-π bonds information |
| `*_Apolar_vdw.csv` | All Apolar VDW bonds information |
| `*_C-H_ON.csv` | All C-H--O/N bonds information |
| `*_C-H_pi.csv` | All C-H--π bonds information |
| `*_Cation_pi.csv` | All Cation-π bonds information |
| `*_Clash.csv` | All Clash bonds information |
| `*_H-bond.csv` | All Hydrogen bond information |
| `*_Halogen_bond.csv` | All Halogen bonds information |
| `*_Lone_pair_pi.csv` | All Lone Pair-π bonds information |
| `*_Metal_Mediated.csv` | All metal Mediated bonds information |
| `*_N-S-O-H_pi.csv` | All O/N/SH-π bonds information |
| `*pi-pi.csv` | All π-π bonds information |
| `*_Polar_vdw.csv` | All Polar VDW bonds information |
| `*_Proximal.csv` | All Proximal bonds information |
| `*_Salt_bridge.csv` | All Salt-bridge bonds information |
| `*_small_summary.csv` | Information about interacting residue count per chain |
| `*_SS_bond.csv` | All SS bonds information |
| `*_summary_table.csv` | Information about count of each type of interaction found |
| `*_Water_Mediated.csv` | All Water-Mediated bonds information |
| `*_ASA_table_chain1.csv` | ASA Table for Chain1 |
| `*_ASA_table_chain2.csv` | ASA Table for Chain2 |
| `*_Rsa_stats.csv` | Information abuot RSA and its stats |

You can find the explanation and information about bonds on the webpage.

There will be some extra files produces during the runtime of code. These files are mainly for processing purposes.

All the output files will be created in the same directory as of the input .pdb file.

---

