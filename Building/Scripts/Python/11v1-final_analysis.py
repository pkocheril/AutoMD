#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:17:29 2025

@author: pkocheril
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# -------------------------
# Helper function to read .xvg files
# -------------------------
def ipt_xvg(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Remove comment/meta lines
    lines = [line.strip() for line in lines
             if not (line.startswith('#') or line.startswith('@') or line.startswith('&'))]

    # Convert to float array
    data = np.array([[float(x) for x in line.split()] for line in lines])
    return data

# -------------------------
# Load atomic charges from .itp
# -------------------------
uv1_idx = [7, 8]  # 1-based indices from MATLAB
uv2_idx = [8, 7]

workingdirectory = os.getcwd()
itplist = [f for f in os.listdir(workingdirectory) if f.endswith('.itp')]

with open(os.path.join(workingdirectory, itplist[0]), 'r') as f:
    A = f.readlines()

line_idx = [None, None]
for i, line in enumerate(A):
    if '[ atoms ]' in line:
        line_idx[0] = i + 2
    elif '[ bonds ]' in line:
        line_idx[1] = i - 2

# Extract atom info
chg = []
atom = []
atom_idx = []

for m in range(line_idx[0], line_idx[1] + 1):
    temp = A[m].strip().split()
    atom_name = temp[5]  # atom name column
    atom.append(atom_name)
    idx = int(''.join([c for c in atom_name if c.isdigit()]))  # extract digits from atom name
    atom_idx.append(idx)
    chg.append(float(temp[6]))  # <-- corrected index for charge

# Wrong index
# chg = []
# atom = []
# atom_idx = []
# for m in range(line_idx[0], line_idx[1] + 1):
#     temp = A[m].split()
#     atom_name = temp[5]
#     atom.append(atom_name)
#     # Extract digits from atom name
#     idx = int(''.join([c for c in atom_name if c.isdigit()]))
#     atom_idx.append(idx)
#     chg.append(float(temp[6]))



chg = np.array(chg)
atom_idx = np.array(atom_idx)

# -------------------------
# Load MD data
# -------------------------
files = glob.glob('*.xvg')
for fname in files:
    if 'xyz.xvg' in fname:
        coord = ipt_xvg(fname)
    elif 'f0.xvg' in fname:
        f0 = ipt_xvg(fname)
    elif 'f1.xvg' in fname:
        f1 = ipt_xvg(fname)

# Remove time column
coord = coord[:, 1:]
f0 = f0[:, 1:]
f1 = f1[:, 1:]

# -------------------------
# Electric field
# -------------------------
ef = f1 - f0

def vecnorm(v, axis=1):
    return np.linalg.norm(v, axis=axis, keepdims=True)

uv1 = coord[:, (uv1_idx[0]-1)*3:(uv1_idx[0]-1)*3+3] - coord[:, (uv1_idx[1]-1)*3:(uv1_idx[1]-1)*3+3]
uv2 = coord[:, (uv2_idx[0]-1)*3:(uv2_idx[0]-1)*3+3] - coord[:, (uv2_idx[1]-1)*3:(uv2_idx[1]-1)*3+3]
uv1 = uv1 / vecnorm(uv1)
uv2 = uv2 / vecnorm(uv2)

# -------------------------
# Create a charge lookup (atom number -> charge)
# -------------------------
charge_map = {atom_idx[i]: chg[i] for i in range(len(atom))}

# -------------------------
# Calculate Fvib1 and Fvib2 with correct charge mapping
# -------------------------
Fvib1 = np.zeros((len(coord), len(atom)))
Fvib2 = np.zeros((len(coord), len(atom)))

for i in range(len(atom)):
    print(f"Atom {atom[i]:<5}  Number={atom_idx[i]:>3}  Charge={charge_map[atom_idx[i]]:>8.4f}")


for n in range(len(atom)):
    eforce = ef[:, (atom_idx[n]-1)*3:(atom_idx[n]-1)*3+3]
    
    # IMPORTANT: use atom number as key in charge_map, not array index
    efield = eforce / charge_map[atom_idx[n]] * 0.1036427
    
    Fvib1[:, n] = np.sum(efield * uv1, axis=1)
    Fvib2[:, n] = np.sum(efield * uv2, axis=1)


# Fvib1 = np.zeros((len(coord), len(atom)))
# Fvib2 = np.zeros((len(coord), len(atom)))

# for n in range(len(atom)):
#     eforce = ef[:, (atom_idx[n]-1)*3:(atom_idx[n]-1)*3+3]
#     efield = eforce / chg[atom_idx[n]-1] * 0.1036427
#     Fvib1[:, n] = np.sum(efield * uv1, axis=1)
#     Fvib2[:, n] = np.sum(efield * uv2, axis=1)

# -------------------------
# Average between atoms
# -------------------------
F1 = np.mean(Fvib1[:, [uv1_idx[0]-1, uv1_idx[1]-1]], axis=1)
F2 = np.mean(Fvib2[:, [uv2_idx[0]-1, uv2_idx[1]-1]], axis=1)

# -------------------------
# Gaussian fit function
# -------------------------
def gaussian_fit(data, histN=128):
    N, edges = np.histogram(data, bins=histN)
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

# -------------------------
# Plot results
# -------------------------
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
plt.show()
