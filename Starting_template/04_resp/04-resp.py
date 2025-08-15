#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  8 21:09:24 2025

@author: pkocheril
"""

from pathlib import Path
import pandas as pd
import numpy as np
import re

# Read the vacuum charges from SP_VAC.chg
with open("Molecule_SP_VAC.chg", "r") as f:
    lines = f.readlines()

# Extract the 5th column from lines with enough elements
vac = []
for line in lines:
    parts = re.split(r"\s+", line.strip())
    if len(parts) >= 5:
        try:
            vac.append(float(parts[4]))
        except ValueError:
            continue
vac = np.array(vac)

# Find all .chg* files
chg_files = sorted(Path(".").glob("*.chg*"))
solvent_names = [f.stem.split("SP_")[1].split(".chg")[0] for f in chg_files]

# Create a DataFrame to hold all RESP data
resp = pd.DataFrame(index=range(len(vac)), columns=solvent_names)

# Process each file
for file_path, name in zip(chg_files, solvent_names):
    with file_path.open() as f:
        lines = f.readlines()

    charges = []
    modified_lines = []

    for line, v in zip(lines, vac):
        parts = re.split(r"\s+", line.strip())
        if len(parts) >= 5:
            try:
                chg = float(parts[4])
                temp = 0.4 * v + 0.6 * chg
                charges.append(temp)
                parts[4] = f"{temp:.6f}"
                modified_lines.append(" ".join(parts))
            except ValueError:
                modified_lines.append(line.strip())
        else:
            modified_lines.append(line.strip())

    resp[name] = charges

    # Save modified RESP file
    with open(f"{name}.resp", "w") as f:
        f.write("\n".join(modified_lines) + "\n")

# Save summary table
resp.to_csv("all.resp", sep=" ", index=False)
