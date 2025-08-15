#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 09:51:05 2025

@author: pkocheril

Generates a .mol2 file from a Gaussian output file.
"""

molname = "Molecule"
OBABEL_PATH = "/central/groups/WeiLab/software/openbabel/bin/obabel"
#OBABEL_PATH = "/usr/local/openbabel/bin/obabel"

import subprocess

def run_cmd(cmd, cwd=None):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=cwd)
    
def convert_to_fchk(inputname):
    run_cmd(f"formchk {inputname}.chk")

def convert_fchk_to_mol2(molname):
    fchk = f"{molname}.fchk"
    mol2 = f"{molname}.mol2"
    run_cmd([OBABEL_PATH, str(fchk), "-O", str(mol2)])
    return mol2

if __name__ == "__main__":
    convert_to_fchk(molname)
    mol2 = convert_fchk_to_mol2(molname)
