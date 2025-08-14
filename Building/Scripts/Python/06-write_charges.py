#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 14:49:56 2025

@author: pkocheril
"""

import pandas as pd
import numpy as np

molname = "Benzonitrile"

# Input files
template_itp = f"{molname}_opt.itp"
resp_file = "all.resp"

# Read RESP data
resp_df = pd.read_csv(resp_file, sep=r"\s+", engine="python")
#resp_df['avg_charge'] = resp_df.mean(axis=1)
columns = ['VAC','water']
resp_df['avg_charge'] = np.average(resp_df[columns], weights=[0.4,0.6], axis=1)

# Read the template .itp file
with open(template_itp, "r") as f:
    lines = f.readlines()

# Locate the [ atoms ] section
start_idx = end_idx = None
for i, line in enumerate(lines):
    if "[ atoms ]" in line:
        start_idx = i + 2  # Skip header lines
    elif "[ bonds ]" in line and start_idx is not None:
        end_idx = i - 2
        break

if start_idx is None or end_idx is None:
    raise ValueError("Could not find [ atoms ] or [ bonds ] sections correctly.")

# Replace charges with 0.4*vacuum+0.6*water

charges = resp_df['avg_charge'].tolist()
mod_lines = lines.copy()

for i, idx in enumerate(range(start_idx, end_idx + 1)):
    parts = mod_lines[idx].split()
    if len(parts) >= 8:
        parts[6] = f"{charges[i]:.6f}"  # 7th column (index 6)
        mod_lines[idx] = " ".join(parts) + "\n"

output_file = f"{molname}-solv.itp"
with open(output_file, "w") as fout:
    fout.writelines(mod_lines)

# names = list(resp_df.columns)

# # Replace charges for each solvent and write to new .itp files
# for name in names:
#     charges = resp_df[name].tolist()
#     mod_lines = lines.copy()

#     for i, idx in enumerate(range(start_idx, end_idx + 1)):
#         parts = mod_lines[idx].split()
#         if len(parts) >= 8:
#             parts[6] = f"{charges[i]:.6f}"  # 7th column (index 6)
#             mod_lines[idx] = " ".join(parts) + "\n"

#     output_file = f"{molname}-{name}.itp"
#     with open(output_file, "w") as fout:
#         fout.writelines(mod_lines)
