#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 12:36:27 2025

@author: pkocheril
"""

import os
import subprocess
from pathlib import Path
from collections import defaultdict
import glob
import argparse
import math

mol_name = "Molecule" # name of ChemDraw file (gets renamed to "Molecule" in bash script)
mol_charge = 0 # 0 by default, will update if needed
mol_spin = 1 # singlets by default
cores = 4 # 4 by default, will use argparse to get nproc
ram = 4 # 4 by default, will use argparse to get nmem (custom script)

# Set paths
if os.path.isdir("/resnick/groups/WeiLab/software/"): # if on the cluster
    OBABEL_PATH = "/resnick/groups/WeiLab/software/openbabel/bin/obabel"
    G16_PATH = "/resnick/groups/WeiLab/software/g16/g16"
    FORMCHK_PATH = "/resnick/groups/WeiLab/software/g16/formchk"
else:
    OBABEL_PATH = "/usr/local/openbabel/bin/obabel"
    G16_PATH = "/Applications/GaussView6/gv/g16"
    FORMCHK_PATH = "/Applications/GaussView6/gv/formchk"



################## Functions ##################

def run_cmd(cmd, cwd=None):
    "Runs a command-line command."
    
    #print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=cwd)

def look_for_files(ext, cwd=None):
    """Searches current or specified directory 
    for files of a specific extension."""
    
    found_files = []
    if cwd:
        found_files = glob.glob(f"{cwd}/*.{ext}")
    else:
        found_files = glob.glob(f"*.{ext}")
    return found_files

def make_gjf(mol_path, job_base_name, charge=0, multiplicity=1, cpus=31, memory=190):
    "Generates a .gjf input for molecule parametrization."
    
    mol_file = Path(mol_path)
    assert mol_file.exists(), f"MOL file not found: {mol_path}"
    
    with mol_file.open() as f:
        mol_lines = f.readlines()

    atom_lines = []
    bonds = defaultdict(list)

    atom_index = 1
    atom_indices = {}  # maps line number to atom index for bonds

    for line in mol_lines:
        parts = line.strip().split()
        if len(parts) >= 16:
            x, y, z = map(float, parts[:3])
            atom = parts[3]
            atom_lines.append(f"{atom:<2} {x:>12.6f} {y:>12.6f} {z:>12.6f}")
            atom_indices[len(atom_lines)] = atom_index
            atom_index += 1
        elif len(parts) == 7 and all(p.isdigit() for p in parts[:3]):
            a1, a2, bond_order = map(int, parts[:3])
            bonds[a1].append((a2, float(bond_order)))

    # Build Gaussian connectivity section
    connectivity_lines = []
    for atom_idx in range(1, len(atom_lines) + 1):
        conn = bonds.get(atom_idx, [])
        if conn:
            line = f"{atom_idx} " + " ".join(f"{i} {b:.1f}" for i, b in conn)
        else:
            line = f"{atom_idx}"
        connectivity_lines.append(line)

    full_input = ""

    jobs = [
        {
            "oldcheck": "",
            "checkpoint": f"{mol_name}_roughopt",
            "method": "HF",
            "basis": "STO-3G",
            "keywords": "opt geom=connectivity", # can't use empirical dispersion here
            "title": f"{mol_name} HF/STO-3G rough geometry optimization",
        },
        {
            "oldcheck": f"{mol_name}_roughopt",
            "checkpoint": f"{mol_name}_optfreq",
            "method": "B3LYP",
            "basis": "6-31G(d,p)",
            "keywords": "opt freq em=GD3BJ geom=allcheck",
            "title": f"{mol_name} B3LYP/6-31G(d,p) opt freq",
        },
        {
            "oldcheck": f"{mol_name}_optfreq",
            "checkpoint": f"{mol_name}_tzvp_opt",
            "method": "B3LYP",
            "basis": "TZVP",
            "keywords": "opt em=GD3BJ geom=allcheck",
            "title": f"{mol_name} B3LYP/TZVP geometry optimization",
        },
        {
            "oldcheck": f"{mol_name}_tzvp_opt",
            "checkpoint": f"{mol_name}_SP_VAC",
            "method": "B3LYP",
            "basis": "def2TZVP",
            "keywords": "em=GD3BJ geom=allcheck",
            "title": f"{mol_name} B3LYP/def2TZVP single point energy in vacuum",
        },
        {
            "oldcheck": f"{mol_name}_tzvp_opt",
            "checkpoint": f"{mol_name}_SP_water",
            "method": "B3LYP",
            "basis": "def2TZVP",
            "keywords": "em=GD3BJ scrf=(solvent=water) geom=allcheck",
            "title": f"{mol_name} B3LYP/def2TZVP single point energy in water",
        },
    ]

    for i, job in enumerate(jobs):
        link = "\n--link1--\n" if i > 0 else ""
        chargemult = f"\n{charge} {multiplicity}\n" if i == 0 else ""
        checks = (
            f"%oldchk={job['oldcheck']}.chk\n%chk={job['checkpoint']}.chk"
            if i > 0 else
            f"%chk={job['checkpoint']}.chk"
        )
        route = f"# {job['method']}/{job['basis']} {job['keywords']}"
        block = (
            f"{link}%nprocshared={cpus}\n"
            f"%mem={memory}GB\n"
            f"{checks}\n"
            f"{route}\n\n"
            f"{job['title']}\n"
            f"{chargemult}"
        )

        if i == 0:
            block += "\n".join(atom_lines) + "\n\n" + "\n".join(connectivity_lines) + "\n"

        block += "\n"
        full_input += block

    input_path = mol_file.parent / f"{job_base_name}_linked.gjf"
    input_path.write_text(full_input)
    return input_path



################## Main ##################

if __name__ == "__main__":
    
    # Make sure directories exist
    input_dir = Path("../00_inputs")
    input_dir.mkdir(exist_ok=True)
    target_dir = Path("../01_QM")
    target_dir.mkdir(exist_ok=True)
    
    # Look for ChemDraw files
    found_files = look_for_files("cdx",input_dir)
    if not found_files:
        found_files = look_for_files("cdxml")
    if found_files:
        input_file = found_files[0]
    else:
        print("ERROR: No ChemDraw files found in the current directory.")
    
    # Convert to SMILES
    smiles_file = f"{target_dir}/{mol_name}.smi"
    run_cmd([OBABEL_PATH, str(input_file), "-O", str(smiles_file)])
    print("Converted to SMILES successfully.")
    
    # Build 3D geometry
    mol3d_file = f"{target_dir}/{mol_name}_3d.mol"
    run_cmd([OBABEL_PATH, str(smiles_file), "-O", str(mol3d_file), "--gen3d"])
    print("Converted to 3D .mol successfully.")
    
    # Coarse optimization
    opt_file = f"{target_dir}/{mol_name}_opt.mol"
    run_cmd([OBABEL_PATH, str(mol3d_file), "-O", str(opt_file), "--minimize"])
    print("Completed coarse geometry optimization successfully.")
    
    # Read total charge via oreport
    run_cmd([OBABEL_PATH, str(input_file), "-oreport", "-O", "oreport.txt"])
    with open("oreport.txt",'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if i == 4:
            charge_line = str(line)
            split_charge = charge_line.split(": ")
            if split_charge[0] == 'TOTAL CHARGE':
                mol_charge = int(split_charge[1])
            else:
                mol_charge = 0
    
    

    ########## Input parsing ##########
    
    # Called as: python 01-qm_prep.py --cores $(nproc) --mem $(nmem)
    parser = argparse.ArgumentParser(description="Determine available computing power.")
    parser.add_argument('--cores', type=int, help='Number of available processing cores.')
    parser.add_argument('--mem', type=float, help='Available memory in GB.')
    args = parser.parse_args()

    # print("Cores:", args.cores)
    # print("Memory:", math.floor(args.mem)," GB")
    if args.cores:
        cores = args.cores
    if args.mem:
        ram = math.floor(args.mem)
    
    
    # Make Gaussian job file
    gjf_path = make_gjf(opt_file, f"{target_dir}/{mol_name}", mol_charge, mol_spin, cores, ram)
    print("Created Gaussian job file successfully.")

    # Run Gaussian job
    run_cmd([G16_PATH, str(gjf_path)])
    print("Gaussian calculation completed successfully.")
    
    # Convert Gaussian outputs
    found_chks = look_for_files("chk")
    for i, file in found_chks:
        file_base = file.strip('.chk')
        run_cmd([FORMCHK_PATH, str(file)])
        run_cmd([OBABEL_PATH, str(f"{file_base}.fchk"), "-O", str(f"{file_base}.mol2")])
        



# def find_cdx_files(cwd=None):
#     cdx_files = []
#     if cwd:
#         os.chdir(cwd)
#     current_directory = os.getcwd()
#     for filename in os.listdir(current_directory):
#         if filename.endswith(".cdx") and os.path.isfile(os.path.join(current_directory, filename)):
#             cdx_files.append(filename)
#     return cdx_files
    # # Look for ChemDraw files
    # found_files = find_cdx_files()
    # if found_files:
    #     input_file = found_files[0]
    # else:
    #     print("ERROR: No .cdx files found in the current directory.")

# 01-build3d_opt.py
# import subprocess
# from pathlib import Path

# molname = "Molecule"

# INPUT_FILE = Path(f"../00_inputs/{molname}.cdx")
# WORKDIR = Path("../01_convert")
# WORKDIR.mkdir(exist_ok=True)
# OPTDIR = Path("../02_geometry")
# OPTDIR.mkdir(exist_ok=True)

# OBABEL_PATH = "/central/groups/WeiLab/software/openbabel/bin/obabel"
# #OBABEL_PATH = "/usr/local/openbabel/bin/obabel"

# def run_cmd(cmd, cwd=None):
#     print(f"Running: {' '.join(cmd)}")
#     subprocess.run(cmd, check=True, cwd=cwd)

# def convert_cdx_to_smiles():
#     smiles = WORKDIR / f"{molname}.smi"
#     run_cmd([OBABEL_PATH, str(INPUT_FILE), "-O", str(smiles)])
#     return smiles

# def generate_3d_structure(smiles_file):
#     mol_3d = WORKDIR / f"{molname}_3d.mol"
#     run_cmd([OBABEL_PATH, str(smiles_file), "-O", str(mol_3d), "--gen3d"])
#     return mol_3d

# def optimize_geometry(input_file):
#     optimized_file = OPTDIR / f"{molname}_opt.mol"
#     run_cmd([OBABEL_PATH, str(input_file), "-O", str(optimized_file), "--minimize"])
#     return optimized_file

# if __name__ == "__main__":
#     print("=== Step 1: CDX to SMILES ===")
#     smiles = convert_cdx_to_smiles()

#     print("=== Step 2: SMILES to 3D structure ===")
#     mol_3d = generate_3d_structure(smiles)

#     print("=== Step 3: Open Babel geometry optimization ===")
#     mol_opt = optimize_geometry(mol_3d)

#     print(f"Completed. Final optimized structure: {mol_opt}")

# 02-gjf_maker.py
# molname = "Molecule"
# molcharge= 0
# molspin = 1
# cores = 12 # 31
# ram = 48 # 190

# from pathlib import Path
# from collections import defaultdict

# def generate_linked_gaussian_input(mol_path, job_base_name, charge=0, multiplicity=1, cpus=31, memory=190):
#     mol_file = Path(mol_path)
#     assert mol_file.exists(), f"MOL file not found: {mol_path}"
    
#     with mol_file.open() as f:
#         mol_lines = f.readlines()

#     atom_lines = []
#     bonds = defaultdict(list)

#     atom_index = 1
#     atom_indices = {}  # maps line number to atom index for bonds

#     for line in mol_lines:
#         parts = line.strip().split()
#         if len(parts) >= 16:
#             x, y, z = map(float, parts[:3])
#             atom = parts[3]
#             atom_lines.append(f"{atom:<2} {x:>12.6f} {y:>12.6f} {z:>12.6f}")
#             atom_indices[len(atom_lines)] = atom_index
#             atom_index += 1
#         elif len(parts) == 7 and all(p.isdigit() for p in parts[:3]):
#             a1, a2, bond_order = map(int, parts[:3])
#             bonds[a1].append((a2, float(bond_order)))

#     # Build Gaussian connectivity section
#     connectivity_lines = []
#     for atom_idx in range(1, len(atom_lines) + 1):
#         conn = bonds.get(atom_idx, [])
#         if conn:
#             line = f"{atom_idx} " + " ".join(f"{i} {b:.1f}" for i, b in conn)
#         else:
#             line = f"{atom_idx}"
#         connectivity_lines.append(line)

#     full_input = ""

#     jobs = [
#         {
#             "oldcheck": "",
#             "checkpoint": f"{molname}_roughopt",
#             "method": "HF",
#             "basis": "STO-3G",
#             "keywords": "opt geom=connectivity", # can't use empirical dispersion here
#             "title": f"{molname} HF/STO-3G rough geometry optimization",
#         },
#         {
#             "oldcheck": f"{molname}_roughopt",
#             "checkpoint": f"{molname}_optfreq",
#             "method": "B3LYP",
#             "basis": "6-31G(d,p)",
#             "keywords": "opt freq em=GD3BJ geom=allcheck",
#             "title": f"{molname} B3LYP/6-31G(d,p) opt freq",
#         },
#         {
#             "oldcheck": f"{molname}_optfreq",
#             "checkpoint": f"{molname}_tzvp_opt",
#             "method": "B3LYP",
#             "basis": "TZVP",
#             "keywords": "opt em=GD3BJ geom=allcheck",
#             "title": f"{molname} B3LYP/TZVP geometry optimization",
#         },
#         {
#             "oldcheck": f"{molname}_tzvp_opt",
#             "checkpoint": f"{molname}_SP_VAC",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ geom=allcheck",
#             "title": f"{molname} B3LYP/def2TZVP single point energy in vacuum",
#         },
#         {
#             "oldcheck": f"{molname}_tzvp_opt",
#             "checkpoint": f"{molname}_SP_water",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ scrf=(solvent=water) geom=allcheck",
#             "title": f"{molname} B3LYP/def2TZVP single point energy in water",
#         },
#     ]

#     for i, job in enumerate(jobs):
#         link = "\n--link1--\n" if i > 0 else ""
#         chargemult = f"\n{charge} {multiplicity}\n" if i == 0 else ""
#         checks = (
#             f"%oldchk={job['oldcheck']}.chk\n%chk={job['checkpoint']}.chk"
#             if i > 0 else
#             f"%chk={job['checkpoint']}.chk"
#         )
#         route = f"# {job['method']}/{job['basis']} {job['keywords']}"
#         block = (
#             f"{link}%nprocshared={cpus}\n"
#             f"%mem={memory}GB\n"
#             f"{checks}\n"
#             f"{route}\n\n"
#             f"{job['title']}\n"
#             f"{chargemult}"
#         )

#         if i == 0:
#             block += "\n".join(atom_lines) + "\n\n" + "\n".join(connectivity_lines) + "\n"

#         block += "\n"
#         full_input += block

#     input_path = mol_file.parent / f"{job_base_name}_linked.gjf"
#     input_path.write_text(full_input)
#     return input_path

# gjf_path = generate_linked_gaussian_input(f"../02_geometry/{molname}_opt.mol", f"../03_gaussian/{molname}", molcharge, molspin, cores, ram)
