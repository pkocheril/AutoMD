#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 09:24:40 2025

@author: pkocheril

Takes an input ChemDraw structure (.cdx or .cdxml) and uses OpenBabel to
convert to SMILES (.smi), build a 3D geometry (_3d.mol), and then do a rough
geometry optimization (_opt.mol)
"""

import subprocess
from pathlib import Path
import glob
import os

inputdir = Path("./00_inputs")


INPUT_FILE = Path("../00_inputs/Benzonitrile.cdx")
WORKDIR = Path("../01_convert")
WORKDIR.mkdir(exist_ok=True)
OPTDIR = Path("../02_geometry")
OPTDIR.mkdir(exist_ok=True)

OBABEL_PATH = "/usr/local/openbabel/bin/obabel"

def find_files_glob(directory, pattern="*"):
    files = glob.glob(os.path.join(directory, pattern))
    return files

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

# Check if inputs provided from previous MD calculation
resp_file = find_files_glob(inputdir, "all.resp") # look for all.resp file

if not resp_file:
    fchk_files = find_files_glob(inputdir, "*.fchk") # look for .fchk files
    
    if not fchk_files:
        cdx_files = find_files_glob(inputdir, "*.cdx") # look for .cdx files
        
        if not cdx_files:
            print("No valid input files found")
            invalidinput = 1
        else:
            for file in cdx_files:
                molname = os.path.basename(file).split(".cdx")[0]

for file in cdx_files:
    print(os.path.basename(file))
    print(os.path.basename(file).split(".cdx")[0])

molname = "Benzonitrile"

if 'molname' in locals():
    print("yes")
    

if __name__ == "__main__":
    print("=== Step 1: CDX to SMILES ===")
    smiles = convert_cdx_to_smiles()

    print("=== Step 2: SMILES to 3D structure ===")
    mol_3d = generate_3d_structure(smiles)

    print("=== Step 3: Open Babel geometry optimization ===")
    mol_opt = optimize_geometry(mol_3d)

    print(f"Completed. Final optimized structure: {mol_opt}")
