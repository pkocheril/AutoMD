#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  8 21:09:24 2025

@author: pkocheril
"""

molname = "Molecule"
molcharge= 0
molspin = 1
cores = 12 # 31
ram = 48 # 190

from pathlib import Path
from collections import defaultdict

def generate_linked_gaussian_input(mol_path, job_base_name, charge=0, multiplicity=1, cpus=31, memory=190):
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
            "checkpoint": f"{molname}_roughopt",
            "method": "HF",
            "basis": "STO-3G",
            "keywords": "opt geom=connectivity", # can't use empirical dispersion here
            "title": f"{molname} HF/STO-3G rough geometry optimization",
        },
        {
            "oldcheck": f"{molname}_roughopt",
            "checkpoint": f"{molname}_optfreq",
            "method": "B3LYP",
            "basis": "6-31G(d,p)",
            "keywords": "opt freq em=GD3BJ geom=allcheck",
            "title": f"{molname} B3LYP/6-31G(d,p) opt freq",
        },
        {
            "oldcheck": f"{molname}_optfreq",
            "checkpoint": f"{molname}_tzvp_opt",
            "method": "B3LYP",
            "basis": "TZVP",
            "keywords": "opt em=GD3BJ geom=allcheck",
            "title": f"{molname} B3LYP/TZVP geometry optimization",
        },
        {
            "oldcheck": f"{molname}_tzvp_opt",
            "checkpoint": f"{molname}_SP_VAC",
            "method": "B3LYP",
            "basis": "def2TZVP",
            "keywords": "em=GD3BJ geom=allcheck",
            "title": f"{molname} B3LYP/def2TZVP single point energy in vacuum",
        },
        {
            "oldcheck": f"{molname}_tzvp_opt",
            "checkpoint": f"{molname}_SP_water",
            "method": "B3LYP",
            "basis": "def2TZVP",
            "keywords": "em=GD3BJ scrf=(solvent=water) geom=allcheck",
            "title": f"{molname} B3LYP/def2TZVP single point energy in water",
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

def create_slurm_submission_script(gaussian_input_path, email="pkocheri@caltech.edu"):
    input_file = Path(gaussian_input_path)
    job_name = input_file.stem
    sh_file = input_file.with_suffix(".sh")

    script = f"""#!/bin/bash

# Submit this script with: sbatch {sh_file.name}

#SBATCH --job-name={job_name}
#SBATCH --time=11:59:59
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-user={email}
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

g16 {input_file.name}
"""

    sh_file.write_text(script)
    return sh_file


gjf_path = generate_linked_gaussian_input(f"../02_geometry/{molname}_opt.mol", f"../03_gaussian/{molname}", molcharge, molspin, cores, ram)
#create_slurm_submission_script(gjf_path)


# from pathlib import Path

# def generate_linked_gaussian_input(mol_path, job_base_name, charge=0, multiplicity=1, cpus=31, memory=190):
#     mol_file = Path(mol_path)
#     assert mol_file.exists(), f"MOL file not found: {mol_path}"
    
#     with mol_file.open() as f:
#         mol_lines = f.readlines()

#     # Extract XYZ from MOL (atom lines only)
#     start_idx = 4
#     atom_lines = []
#     for line in mol_lines[start_idx:]:
#         parts = line.strip().split()
#         if len(parts) > 7:  # Atom line (typically 16 parts)
#             x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
#             atom = parts[3] if parts[3].isalpha() else parts[0]
#             atom_lines.append(f"{atom} {x:.6f} {y:.6f} {z:.6f}")
#         elif len(parts) == 7:  # Bond line — stop parsing
#             break
#     atom_lines.append("\n")

#     if not atom_lines:
#         raise ValueError("No atom coordinates found in MOL file.")

#     # Define job sequence
#     jobs = [
#         {
#             "oldcheck": "",
#             "checkpoint": f"{job_base_name}_roughopt",
#             "method": "HF",
#             "basis": "STO-3G",
#             "keywords": "opt em=GD3BJ",
#             "title": f"{job_base_name} HF/STO-3G rough geometry optimization",
#             "coords": "\n".join(atom_lines),
#         },
#         {
#             "oldcheck": f"{job_base_name}_roughopt",
#             "checkpoint": f"{job_base_name}_optfreq",
#             "method": "B3LYP",
#             "basis": "6-31G(d,p)",
#             "keywords": "opt freq geom=allcheck",
#             "title": f"{job_base_name} B3LYP/6-31G(d,p) opt freq",
#             "coords": "",
#         },
#         {
#             "oldcheck": f"{job_base_name}_optfreq",
#             "checkpoint": f"{job_base_name}_tzvp_opt",
#             "method": "B3LYP",
#             "basis": "TZVP",
#             "keywords": "opt em=GD3BJ geom=allcheck",
#             "title": f"{job_base_name} B3LYP/TZVP geometry optimization",
#             "coords": "",
#         },
#         {
#             "oldcheck": f"{job_base_name}_tzvp_opt",
#             "checkpoint": f"{job_base_name}_SP_VAC",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ geom=allcheck",
#             "title": f"{job_base_name} B3LYP/def2TZVP single point energy in vacuum",
#             "coords": "",
#         },
#         {
#             "oldcheck": f"{job_base_name}_tzvp_opt",
#             "checkpoint": f"{job_base_name}_SP_water",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ scrf=(solvent=water) geom=allcheck",
#             "title": f"{job_base_name} B3LYP/def2TZVP single point energy in water",
#             "coords": "",
#         },
#     ]

#     # Build the full input file
#     full_input = ("")
#     for i, job in enumerate(jobs):
#         link = "\n--link1--\n" if i > 0 else ""
#         checks = f"%oldchk={job['oldcheck']}.chk\n%chk={job['checkpoint']}.chk".strip() if i > 0 else f"%chk={job['checkpoint']}.chk"
#         route = f"# {job['method']}/{job['basis']} {job['keywords']}".strip()
#         block = (
#             f"{link}%nprocshared={cpus}\n"
#             f"%mem={memory}GB\n"
#             f"{checks}\n"
#             f"{route}\n\n"
#             f"{job['title']}\n\n"
#             f"{charge} {multiplicity}\n"
#             f"{job['coords']}\n"
#         )
#         full_input += block

#     input_path = mol_file.parent / f"{job_base_name}_linked.gjf"
#     input_path.write_text(full_input)
#     return input_path

# generate_linked_gaussian_input("water.mol", "water",0,1,4,4)

# from pathlib import Path

# def generate_linked_gaussian_input(mol_path, job_base_name, charge=0, multiplicity=1, cpus=31, memory=190):
#     mol_file = Path(mol_path)
#     assert mol_file.exists(), f"MOL file not found: {mol_path}"
    
#     mol_content = mol_file.read_text().strip()  # Preserve whole MOL contents

#     # Define job sequence
#     jobs = [
#         {
#             "oldcheck": "",
#             "checkpoint": "roughopt",
#             "method": "HF",
#             "basis": "STO-3G",
#             "keywords": "opt geom=connectivity",
#             "title": f"{job_base_name} HF/STO-3G rough geometry optimization",
#             "coords": mol_content,
#         },
#         {
#             "oldcheck": "roughopt",
#             "checkpoint": "optfreq",
#             "method": "B3LYP",
#             "basis": "6-31G(d,p)",
#             "keywords": "opt freq geom=allcheck",
#             "title": f"{job_base_name} B3LYP/6-31G(d,p) opt freq",
#             "coords": "",
#         },
#         {
#             "oldcheck": "optfreq",
#             "checkpoint": "tzvp_opt",
#             "method": "B3LYP",
#             "basis": "TZVP",
#             "keywords": "opt em=GD3BJ geom=allcheck",
#             "title": f"{job_base_name} B3LYP/TZVP geometry optimization",
#             "coords": "",
#         },
#         {
#             "oldcheck": "tzvp_opt",
#             "checkpoint": "SP_VAC",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ geom=allcheck",
#             "title": f"{job_base_name} B3LYP/def2TZVP single point energy in vacuum",
#             "coords": "",
#         },
#         {
#             "oldcheck": "tzvp_opt",
#             "checkpoint": "SP_water",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ scrf=(solvent=water) geom=allcheck",
#             "title": f"{job_base_name} B3LYP/def2TZVP single point energy in water",
#             "coords": "",
#         },
#     ]

#     # Build the full input file
#     full_input = ""
#     for i, job in enumerate(jobs):
#         link = "\n--link1--\n" if i > 0 else "\n"
#         checks = f"%oldchk={job['oldcheck']}.chk\n%chk={job['checkpoint']}.chk".strip() if i > 0 else f"%chk={job['checkpoint']}.chk"
#         route = f"# {job['method']}/{job['basis']} {job['keywords']}".strip()
        
#         coords_block = f"{charge} {multiplicity}\n{job['coords']}\n" if i == 0 else ""
        
#         block = (
#             f"{link}%nprocshared={cpus}\n"
#             f"%mem={memory}GB\n"
#             f"{checks}\n"
#             f"{route}\n\n"
#             f"{job['title']}\n\n"
#             f"{coords_block}\n"
#         )
#         full_input += block

#     input_path = mol_file.parent / f"{job_base_name}_linked.gjf"
#     input_path.write_text(full_input)
#     return input_path

# generate_linked_gaussian_input("water.mol", "water", charge=0, multiplicity=1, cpus=4, memory=4)


# from pathlib import Path

# def generate_linked_gaussian_input(mol_path, job_base_name, charge=0, multiplicity=1, cpus=31, memory=190):
#     mol_file = Path(mol_path)
#     assert mol_file.exists(), f"MOL file not found: {mol_path}"
    
#     with mol_file.open() as f:
#         mol_lines = f.readlines()

#     # Extract XYZ from MOL
#     start_idx = 4
#     atom_lines = []
#     for line in mol_lines[start_idx:]:
#         if line.strip() == '' or len(line.split()) < 4:
#             break
#         parts = line.split()
#         x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
#         atom = parts[3].strip().capitalize() if len(parts) >= 4 and parts[3].strip().isalpha() else "C"
#         atom_lines.append(f"{atom} {x:.6f} {y:.6f} {z:.6f}")

#     if not atom_lines:
#         raise ValueError("No atom coordinates found in MOL file.")

#     # Define job sequence
#     jobs = [
#         {
#             "oldcheck": "",
#             "checkpoint": f"{job_base_name}_roughopt",
#             "method": "HF",
#             "basis": "STO-3G",
#             "keywords": "opt geom=connectivity",
#             "title": f"{job_base_name} HF/STO-3G rough geometry optimization",
#             "coords": "\n".join(atom_lines),
#         },
#         {
#             "oldcheck": f"{job_base_name}_roughopt",
#             "checkpoint": f"{job_base_name}_optfreq",
#             "method": "B3LYP",
#             "basis": "6-31G(d,p)",
#             "keywords": "opt freq geom=allcheck",
#             "title": f"{job_base_name} B3LYP/6-31G(d,p) opt freq",
#             "coords": "",
#         },
#         {
#             "oldcheck": f"{job_base_name}_optfreq",
#             "checkpoint": f"{job_base_name}_tzvp_opt",
#             "method": "B3LYP",
#             "basis": "TZVP",
#             "keywords": "opt em=GD3BJ geom=allcheck",
#             "title": f"{job_base_name} B3LYP/TZVP geometry optimization",
#             "coords": "",
#         },
#         {
#             "oldcheck": f"{job_base_name}_tzvp_opt",
#             "checkpoint": f"{job_base_name}_SP_VAC",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ geom=allcheck",
#             "title": f"{job_base_name} B3LYP/def2TZVP single point energy in vacuum",
#             "coords": "",
#         },
#         {
#             "oldcheck": f"{job_base_name}_tzvp_opt",
#             "checkpoint": f"{job_base_name}_SP_water",
#             "method": "B3LYP",
#             "basis": "def2TZVP",
#             "keywords": "em=GD3BJ scrf=(solvent=water) geom=allcheck",
#             "title": f"{job_base_name} B3LYP/def2TZVP single point energy in water",
#             "coords": "",
#         },
#     ]

#     # Build the full input file
#     full_input = ("")
#     for i, job in enumerate(jobs):
#         link = "\n--link1--\n" if i > 0 else "\n"
#         checks = f"%oldchk={job['oldcheck']}.chk\n%chk={job['checkpoint']}.chk".strip() if i > 0 else f"%chk={job['checkpoint']}.chk"
#         route = f"# {job['method']}/{job['basis']} {job['keywords']}".strip()
#         block = (
#             f"{link}%nprocshared={cpus}\n"
#             f"%mem={memory}GB\n"
#             f"{checks}\n"
#             f"{route}\n\n"
#             f"{job['title']}\n\n"
#             f"{charge} {multiplicity}\n"
#             f"{job['coords']}\n"
#             f"\n"
#         )
#         full_input += block

#     input_path = mol_file.parent / f"{job_base_name}_linked.gjf"
#     input_path.write_text(full_input)
#     return input_path

# generate_linked_gaussian_input("water.mol", "water",0,1,4,4)

