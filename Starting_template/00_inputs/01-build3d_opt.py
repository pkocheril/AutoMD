#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  8 21:09:24 2025

@author: pkocheril
"""

import subprocess
from pathlib import Path

molname = "Molecule"

INPUT_FILE = Path(f"../00_inputs/{molname}.cdx")
WORKDIR = Path("../01_convert")
WORKDIR.mkdir(exist_ok=True)
OPTDIR = Path("../02_geometry")
OPTDIR.mkdir(exist_ok=True)

OBABEL_PATH = "/central/groups/WeiLab/software/openbabel/bin/obabel"
#OBABEL_PATH = "/usr/local/openbabel/bin/obabel"

def run_cmd(cmd, cwd=None):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=cwd)

def convert_cdx_to_smiles():
    smiles = WORKDIR / f"{molname}.smi"
    run_cmd([OBABEL_PATH, str(INPUT_FILE), "-O", str(smiles)])
    return smiles

def generate_3d_structure(smiles_file):
    mol_3d = WORKDIR / f"{molname}_3d.mol"
    run_cmd([OBABEL_PATH, str(smiles_file), "-O", str(mol_3d), "--gen3d"])
    return mol_3d

def optimize_geometry(input_file):
    optimized_file = OPTDIR / f"{molname}_opt.mol"
    run_cmd([OBABEL_PATH, str(input_file), "-O", str(optimized_file), "--minimize"])
    return optimized_file

if __name__ == "__main__":
    print("=== Step 1: CDX to SMILES ===")
    smiles = convert_cdx_to_smiles()

    print("=== Step 2: SMILES to 3D structure ===")
    mol_3d = generate_3d_structure(smiles)

    print("=== Step 3: Open Babel geometry optimization ===")
    mol_opt = optimize_geometry(mol_3d)

    print(f"Completed. Final optimized structure: {mol_opt}")


# #!/usr/bin/env python3
# import subprocess
# from pathlib import Path

# INPUT_FILE = Path("inputs/Benzonitrile.cdx")
# WORKDIR = Path("stages/00_preparation")
# WORKDIR.mkdir(exist_ok=True)

# OBABEL_PATH = "/usr/local/openbabel/bin/obabel"
# OBMIN_PATH = "/usr/local/openbabel/bin/obminimize"

# def run_cmd(cmd, cwd=None):
#     print(f"Running: {' '.join(cmd)}")
#     subprocess.run(cmd, check=True, cwd=cwd)

# def convert_cdx_to_smiles():
#     smiles = WORKDIR / "molecule.smi"
#     run_cmd([OBABEL_PATH, str(INPUT_FILE), "-O", str(smiles)])
#     return smiles

# def generate_3d_structure(smiles_file):
#     mol_3d = WORKDIR / "molecule_3d.mol"
#     run_cmd([OBABEL_PATH, str(smiles_file), "-O", str(mol_3d), "--gen3d"])
#     return mol_3d

# def geometry_optimize(mol_3d):
#     opt_mol = WORKDIR / "molecule_optimized.mol"
#     run_cmd([OBABEL_PATH, str(mol_3d), "-O", str(opt_mol), "--minimize"])

# if __name__ == "__main__":
#     print("=== Step 1: CDX to SMILES ===")
#     smiles = convert_cdx_to_smiles()

#     print("=== Step 2: SMILES to 3D structure ===")
#     mol_3d = generate_3d_structure(smiles)
    
#     print("=== Step 3: Geometry optimization ===")
#     opt_mol = geometry_optimize(mol_3d)

#     print(f"Completed. Final output: {mol_3d}")
