#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:36:14 2025

@author: pkocheril
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# ================================
# Constants
# ================================
MV_CM_CONVERSION = 0.1036427
HIST_BINS = 128

# ================================
# File parsing functions
# ================================
def read_itp(filename):
    """
    Parse an .itp file to extract atom names, atom numbers, and charges.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find section bounds
    line_idx = [None, None]
    for i, line in enumerate(lines):
        if '[ atoms ]' in line:
            line_idx[0] = i + 2
        elif '[ bonds ]' in line:
            line_idx[1] = i - 2

    atom = []
    atom_idx = []
    chg = []
    for m in range(line_idx[0], line_idx[1] + 1):
        temp = lines[m].strip().split()
        atom_name = temp[5]
        atom.append(atom_name)
        idx = int(''.join([c for c in atom_name if c.isdigit()]))
        atom_idx.append(idx)
        chg.append(float(temp[6]))  # fixed charge index

    atom_idx = np.array(atom_idx, dtype=int)
    chg = np.array(chg, dtype=float)
    return atom, atom_idx, chg


def ipt_xvg(filename):
    """
    Load an .xvg file, skipping comments/meta lines.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    lines = [line.strip() for line in lines
             if not (line.startswith('#') or line.startswith('@') or line.startswith('&'))]

    data = np.array([[float(x) for x in line.split()] for line in lines])
    return data


def load_xvg_files():
    """
    Load coord, f0, f1 data from .xvg files in the current directory.
    """
    files = glob.glob('*.xvg')
    coord = f0 = f1 = None

    for fname in files:
        if 'xyz.xvg' in fname:
            coord = ipt_xvg(fname)
        elif 'f0.xvg' in fname:
            f0 = ipt_xvg(fname)
        elif 'f1.xvg' in fname:
            f1 = ipt_xvg(fname)

    if coord is None or f0 is None or f1 is None:
        raise FileNotFoundError("Missing one or more of xyz.xvg, f0.xvg, or f1.xvg")

    # Remove time column
    coord = coord[:, 1:]
    f0 = f0[:, 1:]
    f1 = f1[:, 1:]
    return coord, f0, f1

# ================================
# Computation functions
# ================================
def compute_fields(coord, f0, f1, atom_idx, chg, uv1_idx, uv2_idx):
    """
    Compute projected electric fields along two molecular axes.
    """
    # Electric field = (f1 - f0) / charge
    ef = f1 - f0

    # Build unit vectors
    uv1 = coord[:, (uv1_idx[0]-1)*3:(uv1_idx[0]-1)*3+3] - coord[:, (uv1_idx[1]-1)*3:(uv1_idx[1]-1)*3+3]
    uv2 = coord[:, (uv2_idx[0]-1)*3:(uv2_idx[0]-1)*3+3] - coord[:, (uv2_idx[1]-1)*3:(uv2_idx[1]-1)*3+3]
    uv1 /= np.linalg.norm(uv1, axis=1, keepdims=True)
    uv2 /= np.linalg.norm(uv2, axis=1, keepdims=True)

    # Create charge lookup
    charge_map = {atom_idx[i]: chg[i] for i in range(len(atom_idx))}

    # Vectorized computation: (frames, atoms, 3)
    efields = np.stack([
        ef[:, (atom_idx[n]-1)*3:(atom_idx[n]-1)*3+3] / charge_map[atom_idx[n]] * MV_CM_CONVERSION
        for n in range(len(atom_idx))
    ], axis=1)

    Fvib1 = np.einsum('nij,nj->ni', efields, uv1)
    Fvib2 = np.einsum('nij,nj->ni', efields, uv2)

    # Average between given atoms
    F1 = np.mean(Fvib1[:, [uv1_idx[0]-1, uv1_idx[1]-1]], axis=1)
    F2 = np.mean(Fvib2[:, [uv2_idx[0]-1, uv2_idx[1]-1]], axis=1)
    return F1, F2


def gaussian_fit(data, hist_bins=HIST_BINS):
    """
    Fit a Gaussian to histogrammed data.
    """
    N, edges = np.histogram(data, bins=hist_bins)
    edges = edges[1:] - (edges[1] - edges[0]) / 2

    def fun(r):
        return r[0] * np.exp(-((edges - r[1]) / r[2])**2) - N

    r0 = [np.median(N), 0, 10]
    lb = [0, -np.inf, 0]

    res = least_squares(fun, r0, bounds=(lb, [np.inf, np.inf, np.inf]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=1_000_000)
    fitval = res.x
    fcurve = fitval[0] * np.exp(-((edges - fitval[1]) / fitval[2])**2)
    return edges, N, fcurve, fitval

# ================================
# Main analysis routine
# ================================
def main():
    # Atom pairs (MATLAB 1-based indexing)
    uv1_idx = [7, 8]
    uv2_idx = [8, 7]

    # Read input files
    itp_files = [f for f in os.listdir() if f.endswith('.itp')]
    if not itp_files:
        raise FileNotFoundError("No .itp file found in current directory")
    atom, atom_idx, chg = read_itp(itp_files[0])
    coord, f0, f1 = load_xvg_files()

    # Compute electric fields
    F1, F2 = compute_fields(coord, f0, f1, atom_idx, chg, uv1_idx, uv2_idx)

    # Fit and plot
    edges1, N1, fcurve1, fitval1 = gaussian_fit(F1)
    edges2, N2, fcurve2, fitval2 = gaussian_fit(F2)

    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.plot(edges1, N1, 'r.')
    plt.plot(edges1, fcurve1, 'b')
    plt.title(f"Between {atom[uv1_idx[0]-1]} and {atom[uv1_idx[1]-1]}\n"
              f"{fitval1[1]:.2f} ± {fitval1[2]/np.sqrt(2):.2f} MV/cm")
    plt.xlabel('Electric field (MV/cm)')
    plt.ylabel('Counts')

    plt.subplot(1, 2, 2)
    plt.plot(edges2, N2, 'r.')
    plt.plot(edges2, fcurve2, 'b')
    plt.title(f"Between {atom[uv2_idx[0]-1]} and {atom[uv2_idx[1]-1]}\n"
              f"{fitval2[1]:.2f} ± {fitval2[2]/np.sqrt(2):.2f} MV/cm")
    plt.xlabel('Electric field (MV/cm)')
    plt.ylabel('Counts')

    plt.tight_layout()
    plt.savefig("electric_field_histograms.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
