#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 17:13:26 2025

@author: pkocheril
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def ipt_xvg(filename):
    """
    Reads an .xvg file, removing lines starting with #, @, or &,
    and returns a NumPy array of numeric data.
    """
    if not os.path.exists(filename):
        print(f"Warning: File not found: {filename}")
        return None

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Remove metadata/comment lines
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

# Find all .xvg files in current directory
xvg_files = sorted(glob.glob("*.xvg"))

if not xvg_files:
    print("No .xvg files found in current directory. Exiting.")
else:
    print(f"Found {len(xvg_files)} .xvg files: {xvg_files}")

    # Assign line styles in a repeating cycle
    styles = ['.-', '.', '-', '--']
    datasets = []
    for i, file in enumerate(xvg_files):
        data = ipt_xvg(file)
        if data is not None:
            title = os.path.splitext(os.path.basename(file))[0]
            style = styles[i % len(styles)]
            datasets.append((data, title, style))

    if not datasets:
        print("No valid data to plot. Exiting.")
    else:
        fig, axes = plt.subplots(1, len(datasets), figsize=(5 * len(datasets), 4))

        # Ensure axes is iterable
        if len(datasets) == 1:
            axes = [axes]

        for ax, (data, title, style) in zip(axes, datasets):
            avg_last_half = np.mean(data[len(data)//2:, 1])
            ax.plot(data[:, 0], data[:, 1], style)
            ax.legend([f'avg(last half) = {avg_last_half:.3f}'])
            ax.set_title(title)

        plt.tight_layout()

        # Ensure output directory exists
        save_path = "../08_final_results/solvent_equilibration.png"
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        # Save and show
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to: {save_path}")

        plt.show()
