#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 17:13:26 2025

@author: pkocheril
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys

def ipt_xvg(filename):
    """Reads an .xvg file, removing lines starting with #, @, or &."""
    if not os.path.exists(filename):
        print(f"Warning: File not found: {filename}")
        return None

    with open(filename, 'r') as f:
        lines = f.readlines()

    data_lines = [
        line.strip() for line in lines
        if not line.startswith(('#', '@', '&')) and line.strip()
    ]

    try:
        data = np.array([list(map(float, line.split())) for line in data_lines])
    except ValueError as e:
        print(f"Error reading numeric data from {filename}: {e}")
        return None

    return data

# Find all .xvg files
xvg_files = sorted(glob.glob("*.xvg"))

if not xvg_files:
    print("Error: No .xvg files found.")
    sys.exit(1)  # Exit code 1 = Missing input files

styles = ['.-', '.', '-', '--']
datasets = []
for i, file in enumerate(xvg_files):
    data = ipt_xvg(file)
    if data is not None:
        title = os.path.splitext(os.path.basename(file))[0].capitalize()
        style = styles[i % len(styles)]
        datasets.append((data, title, style))

if not datasets:
    print("Error: All files failed to load valid data.")
    sys.exit(2)  # Exit code 2 = Files found but unreadable

fig, axes = plt.subplots(1, len(datasets), figsize=(5 * len(datasets), 4))
if len(datasets) == 1:
    axes = [axes]

for ax, (data, title, style) in zip(axes, datasets):
    avg_last_half = np.mean(data[len(data)//2:, 1])
    ax.plot(data[:, 0], data[:, 1], style)
    ax.legend([f'avg(last half) = {avg_last_half:.3f}'])
    ax.set_title(title)

plt.tight_layout()

save_path = "../06_final_results/solvent_equilibration.png"
os.makedirs(os.path.dirname(save_path), exist_ok=True)
plt.savefig(save_path, dpi=300)
plt.close(fig)

print(f"Plot saved to: {save_path}")
sys.exit(0)  # Exit code 0 = Success
