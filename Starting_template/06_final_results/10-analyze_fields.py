#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:57:44 2025

@author: pkocheril
"""

import os
import glob
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# ================================
# Constants & thresholds
# ================================
MV_CM_CONVERSION = 0.1036427
HIST_BINS = 128

# Bond-length (Å) heuristics for auto-detection (on time-averaged structure)
NITRILE_MIN, NITRILE_MAX = 1.10, 1.22   # C≡N
CARBONYL_MIN, CARBONYL_MAX = 1.18, 1.32 # C=O

# ================================
# Utilities
# ================================
def element_from_name(name: str) -> str:
    """Best-effort element guess from atom name (first letter)."""
    return ''.join([c for c in name if c.isalpha()])[:1].upper() or '?'

def avg_positions_from_coord(coord: np.ndarray) -> np.ndarray:
    """Return (n_atoms, 3) average positions from coord matrix (frames, 3N)."""
    n_atoms = coord.shape[1] // 3
    xyz = coord.reshape(coord.shape[0], n_atoms, 3)
    return xyz.mean(axis=0)

def distances_for_bonds(avg_pos: np.ndarray, bonds: list[tuple[int,int]]) -> dict[tuple[int,int], float]:
    d = {}
    for a, b in bonds:
        va = avg_pos[a-1]
        vb = avg_pos[b-1]
        d[(a, b)] = float(np.linalg.norm(va - vb))
    return d

# ================================
# File parsing
# ================================
def read_itp(filename):
    """
    Parse an .itp file: atoms (number, name, charge) and bonds (pairs of atom numbers).
    Returns:
      atom_names (list[str]),
      atom_numbers (np.ndarray[int]),
      charges (np.ndarray[float]),
      bonds (list[tuple[int,int]])
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Section indices
    atoms_start = atoms_end = bonds_start = bonds_end = None
    for i, line in enumerate(lines):
        if '[ atoms ]' in line:
            atoms_start = i + 2  # typical spacing: header + one blank/comment
        if '[ bonds ]' in line:
            bonds_start = i + 2
            if atoms_end is None:
                atoms_end = i - 2

    # Find end of bonds (next section or EOF)
    if bonds_start is not None:
        for j in range(bonds_start, len(lines)):
            if '[' in lines[j] and ']' in lines[j] and j > bonds_start:
                bonds_end = j - 2
                break
        if bonds_end is None:
            bonds_end = len(lines) - 1

    if atoms_start is None or atoms_end is None:
        raise ValueError("Could not locate a valid [ atoms ] section in the .itp")

    atom_names, atom_numbers, charges = [], [], []
    for m in range(atoms_start, atoms_end + 1):
        row = lines[m].strip()
        if not row or row.startswith(('#',';','@','&')):
            continue
        temp = row.split()

        # GROMACS [ atoms ] canonical columns (most forcefields):
        # nr  type  resnr  resid  atom  cgnr  charge  mass  (…)
        # We’ll try robust parsing:
        try:
            nr = int(temp[0])                   # atom number in topology
            name = temp[4]                      # atom name
            charge = float(temp[6])             # charge (you fixed this earlier)
        except (ValueError, IndexError):
            # Fallback to previous behavior if layout differs:
            name = temp[5]
            nr = int(''.join([c for c in name if c.isdigit()]) or '0')
            charge = float(temp[6])
        atom_numbers.append(nr)
        atom_names.append(name)
        charges.append(charge)

    atom_numbers = np.asarray(atom_numbers, dtype=int)
    charges = np.asarray(charges, dtype=float)

    bonds = []
    if bonds_start is not None and bonds_end is not None:
        for m in range(bonds_start, bonds_end + 1):
            row = lines[m].strip()
            if not row or row.startswith(('#',';','@','&','[')):
                continue
            parts = row.split()
            # First two columns are the bonded atom numbers in most .itp formats
            try:
                a = int(parts[0]); b = int(parts[1])
                bonds.append((a, b) if a < b else (b, a))
            except (ValueError, IndexError):
                continue

    return atom_names, atom_numbers, charges, bonds

def ipt_xvg(filename):
    """Load an .xvg file, skipping comments/meta lines."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines
             if not (line.startswith('#') or line.startswith('@') or line.startswith('&') or not line.strip())]
    data = np.array([[float(x) for x in line.split()] for line in lines])
    return data

def load_xvg_files():
    """Load coord, f0, f1 (strip time col)."""
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
        raise FileNotFoundError("Missing one or more of xyz.xvg, f0.xvg, f1.xvg")
    return coord[:,1:], f0[:,1:], f1[:,1:]

# ================================
# Field computation & fitting
# ================================
def compute_fields(coord, f0, f1, atom_idx, chg, pair):
    """Projected fields along the axis defined by `pair` (1-based atom numbers)."""
    ef = f1 - f0
    a1, a2 = pair
    uv = coord[:, (a1-1)*3:(a1-1)*3+3] - coord[:, (a2-1)*3:(a2-1)*3+3]
    uv /= np.linalg.norm(uv, axis=1, keepdims=True)

    charge_map = {atom_idx[i]: chg[i] for i in range(len(atom_idx))}
    efields = np.stack([
        ef[:, (atom_idx[n]-1)*3:(atom_idx[n]-1)*3+3] / charge_map[atom_idx[n]] * MV_CM_CONVERSION
        for n in range(len(atom_idx))
    ], axis=1)  # (frames, atoms, 3)

    Fvib = np.einsum('nij,nj->ni', efields, uv)   # (frames, atoms)
    Fmean = np.mean(Fvib[:, [a1-1, a2-1]], axis=1)
    return Fmean

def gaussian_fit(data, hist_bins=HIST_BINS):
    """Gaussian fit to histogrammed data."""
    N, edges = np.histogram(data, bins=hist_bins)
    edges = edges[1:] - (edges[1] - edges[0]) / 2
    def fun(r): return r[0]*np.exp(-((edges - r[1]) / r[2])**2) - N
    r0 = [np.median(N), 0, 10]; lb = [0, -np.inf, 0]
    res = least_squares(fun, r0, bounds=(lb, [np.inf, np.inf, np.inf]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=1_000_000)
    fitval = res.x
    fcurve = fitval[0]*np.exp(-((edges - fitval[1]) / fitval[2])**2)
    return edges, N, fcurve, fitval

# ================================
# Auto-detection of nitriles/carbonyls
# ================================
# Add failsafe in case no pairs are found
def guess_pairs(atom_names, atom_numbers, bonds, coord,
                nitrile_range=(NITRILE_MIN, NITRILE_MAX),
                carbonyl_range=(CARBONYL_MIN, CARBONYL_MAX)):
    """
    Guess nitrile (C≡N) and carbonyl (C=O) pairs using connectivity + average bond lengths.
    If none meet the distance criteria, fallback to closest C-N and C-O bonds.
    """
    avg_pos = avg_positions_from_coord(coord)
    bond_len = distances_for_bonds(avg_pos, bonds)
    deg = {n: 0 for n in atom_numbers}
    for a, b in bonds:
        deg[a] += 1; deg[b] += 1
    elem = {atom_numbers[i]: element_from_name(atom_names[i]) for i in range(len(atom_numbers))}

    nitriles, carbonyls = [], []

    for (a, b), d in bond_len.items():
        ea, eb = elem.get(a, '?'), elem.get(b, '?')
        pair = (a, b) if a < b else (b, a)
        # Nitrile detection
        if {ea, eb} == {'C', 'N'} and nitrile_range[0] <= d <= nitrile_range[1]:
            if deg[a] == 1 or deg[b] == 1:
                nitriles.append(pair)
        # Carbonyl detection
        if {ea, eb} == {'C', 'O'} and carbonyl_range[0] <= d <= carbonyl_range[1]:
            carbonyls.append(pair)

    # Fallbacks if no pairs detected
    if not nitriles:
        # pick C-N bond with shortest distance
        cn_bonds = [(pair, d) for pair, d in bond_len.items() if set([elem[pair[0]], elem[pair[1]]]) == {'C','N'}]
        if cn_bonds:
            closest = min(cn_bonds, key=lambda x: x[1])
            nitriles.append(closest[0])
    if not carbonyls:
        # pick C-O bond with shortest distance
        co_bonds = [(pair, d) for pair, d in bond_len.items() if set([elem[pair[0]], elem[pair[1]]]) == {'C','O'}]
        if co_bonds:
            closest = min(co_bonds, key=lambda x: x[1])
            carbonyls.append(closest[0])

    # Deduplicate
    def dedup(seq): 
        seen = set(); out = []
        for x in seq:
            if x not in seen:
                out.append(x); seen.add(x)
        return out

    return {
        'nitrile': dedup(nitriles),
        'carbonyl': dedup(carbonyls),
    }


# ================================
# Main
# ================================
def main():
    parser = argparse.ArgumentParser(description="Electric field analysis with auto-detection of nitriles/carbonyls.")
    parser.add_argument("--pairs", nargs="*", default=None,
                        help="Atom pairs to analyze, e.g., 7,8 8,7 10,15. If omitted, auto-detect pairs.")
    parser.add_argument("--hist-bins", type=int, default=HIST_BINS, help="Histogram bins (default 128).")
    parser.add_argument("--save-prefix", default="electric_field", help="Prefix for output files.")
    parser.add_argument("--nitrile-range", default=f"{NITRILE_MIN},{NITRILE_MAX}",
                        help="Å range for nitrile detection, e.g., 1.10,1.22")
    parser.add_argument("--carbonyl-range", default=f"{CARBONYL_MIN},{CARBONYL_MAX}",
                        help="Å range for carbonyl detection, e.g., 1.18,1.32")
    args = parser.parse_args()

    nitrile_range = tuple(map(float, args.nitrile_range.split(',')))
    carbonyl_range = tuple(map(float, args.carbonyl_range.split(',')))

    # Load inputs
    itp_files = [f for f in os.listdir() if f.endswith('.itp')]
    if not itp_files:
        raise FileNotFoundError("No .itp file found in current directory.")
    atom_names, atom_idx, chg, bonds = read_itp(itp_files[0])
    coord, f0, f1 = load_xvg_files()

    # Decide pairs
    if args.pairs:
        pairs = [tuple(map(int, p.split(","))) for p in args.pairs]
    else:
        guessed = guess_pairs(atom_names, atom_idx, bonds, coord,
                              nitrile_range=nitrile_range, carbonyl_range=carbonyl_range)
        # Prioritize nitriles, then carbonyls
        pairs = guessed['nitrile'] + guessed['carbonyl']
        if not pairs:
            raise RuntimeError("Auto-detection found no nitrile or carbonyl pairs. "
                               "Try adjusting --nitrile-range/--carbonyl-range or pass --pairs manually.")

    # Prepare plotting
    plt.figure(figsize=(5 * len(pairs), 4))
    csv_rows = []

    for i, pair in enumerate(pairs, 1):
        F = compute_fields(coord, f0, f1, atom_idx, chg, pair)
        edges, N, fcurve, fitval = gaussian_fit(F, hist_bins=args.hist_bins)

        # Safe atom labels (if indices exceed table length, print indices)
        def label(atom_no):
            return atom_names[atom_no-1] if 1 <= atom_no <= len(atom_names) else f"Atom{atom_no}"
        label_pair = f"{label(pair[0])}-{label(pair[1])}"

        plt.subplot(1, len(pairs), i)
        plt.plot(edges, N, 'r.')
        plt.plot(edges, fcurve, 'b')
        plt.title(f"{label_pair}\n{fitval[1]:.2f} ± {fitval[2]/np.sqrt(2):.2f} MV/cm")
        plt.xlabel('Electric field (MV/cm)')
        plt.ylabel('Counts')

        csv_rows.append([
            label_pair, pair[0], pair[1],
            fitval[0], fitval[1], fitval[2], fitval[2]/np.sqrt(2)
        ])

    plt.tight_layout()
    fig_path = f"{args.save_prefix}_histograms.png"
    plt.savefig(fig_path, dpi=300)
    plt.show()

    # Save CSV
    csv_path = f"{args.save_prefix}_fit_results.csv"
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Pair(Label)", "Atom1", "Atom2",
            "Amplitude", "Center (MV/cm)", "Sigma (MV/cm)", "Sigma/sqrt(2) (MV/cm)"
        ])
        writer.writerows(csv_rows)

    print(f"Processed {len(pairs)} pairs.")
    print(f"Saved plot: {fig_path}")
    print(f"Saved CSV:  {csv_path}")

if __name__ == "__main__":
    main()
