#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 09:25:48 2025

@author: pkocheril
"""

import os
import glob
# import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import pandas as pd

# ================================
# Constants & thresholds
# ================================
MV_CM_CONVERSION = -0.1036427 # negative to get the sign convention
HIST_BINS = 128

# Heuristics in Å (fallback geometry)
NITRILE_MIN_A, NITRILE_MAX_A = 1.10, 1.22   # C≡N
CARBONYL_MIN_A, CARBONYL_MAX_A = 1.18, 1.32 # C=O

# Convert to nm for .itp [ bonds ] (r0 is in nm)
NITRILE_MIN_NM, NITRILE_MAX_NM = NITRILE_MIN_A * 0.1, NITRILE_MAX_A * 0.1
CARBONYL_MIN_NM, CARBONYL_MAX_NM = CARBONYL_MIN_A * 0.1, CARBONYL_MAX_A * 0.1

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
    """Return distances (Å) for each bond pair (1-based indices)."""
    d = {}
    for a, b in bonds:
        va = avg_pos[a-1]
        vb = avg_pos[b-1]
        d[(a, b)] = float(np.linalg.norm(va - vb))
    return d

def dedup_pairs(seq):
    seen = set(); out = []
    for x in seq:
        if x not in seen:
            out.append(x); seen.add(x)
    return out

# ================================
# File parsing
# ================================
def read_itp(filename):
    """
    Parse an .itp file:
      - [ atoms ] -> atom_numbers, atom_names, charges
      - [ bonds ] -> bonds (pairs) and bonds_meta with r0 (nm) + label from comment
    Returns:
      atom_names (list[str]),
      atom_numbers (np.ndarray[int]),
      charges (np.ndarray[float]),
      bonds (list[tuple[int,int]]),
      bonds_meta (list[dict] with keys: i, j, r0_nm, label)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Locate sections
    atoms_start = atoms_end = bonds_start = bonds_end = None
    for i, line in enumerate(lines):
        if '[ atoms' in line:
            atoms_start = i + 1
        if '[ bonds' in line:
            bonds_start = i + 1
            if atoms_end is None and atoms_start is not None:
                atoms_end = i - 1

    # Find end of bonds (next section or EOF)
    if bonds_start is not None:
        bonds_end = len(lines) - 1
        for j in range(bonds_start, len(lines)):
            if j == bonds_start: 
                continue
            if '[' in lines[j] and ']' in lines[j]:
                bonds_end = j - 1
                break

    if atoms_start is None:
        raise ValueError("Could not locate a valid [ atoms ] section in the .itp")
    if atoms_end is None:
        # until next section or EOF
        atoms_end = len(lines) - 1
        for j in range(atoms_start, len(lines)):
            if j == atoms_start:
                continue
            if '[ ' in lines[j] and ']' in lines[j]:
                atoms_end = j - 1
                break

    # --- Parse [ atoms ] ---
    atom_names, atom_numbers, charges = [], [], []
    for m in range(atoms_start, atoms_end + 1):
        row = lines[m].strip()
        if not row or row.startswith(('#',';','@','&','[')):
            continue
        temp = row.split()
        # Expected (most forcefields):
        # nr  type  resnr  resname  atom  cgnr  charge  mass ...
        try:
            nr = int(temp[0])
            name = temp[4]
            charge = float(temp[6])
        except (ValueError, IndexError):
            # Fallback: try alternate offsets (rare)
            # (Keep this conservative to avoid mis-parsing)
            continue
        atom_numbers.append(nr)
        atom_names.append(name)
        charges.append(charge)

    atom_numbers = np.asarray(atom_numbers, dtype=int)
    charges = np.asarray(charges, dtype=float)

    # --- Parse [ bonds ] ---
    bonds = []
    bonds_meta = []
    if bonds_start is not None and bonds_end is not None:
        for m in range(bonds_start, bonds_end + 1):
            line = lines[m].strip()
            if not line or line.startswith(('#',';','@','&','[')):
                continue

            # Split off comment (everything after ';')
            if ';' in line:
                data_part, comment_part = line.split(';', 1)
                comment = comment_part.strip()
            else:
                data_part, comment = line, ''

            parts = data_part.split()
            if len(parts) < 2:
                continue

            # First two columns are atom indices
            try:
                i_atom = int(parts[0]); j_atom = int(parts[1])
            except ValueError:
                continue
            pair = (i_atom, j_atom) if i_atom < j_atom else (j_atom, i_atom)
            bonds.append(pair)

            # 4th column r0 (nm) if present
            r0_nm = None
            if len(parts) >= 4:
                try:
                    r0_nm = float(parts[3])
                except ValueError:
                    r0_nm = None

            # Extract a simple bond label like "C-N", "N-C", "C-O", "O-C" from the comment
            label = None
            if comment:
                c = comment.upper()
                # be permissive about punctuation around the token
                if 'C-N' in c or 'N-C' in c:
                    label = 'C-N'
                elif 'C-O' in c or 'O-C' in c:
                    label = 'C-O'
                # (extend here if you ever want C=C, etc.)

            bonds_meta.append({'i': pair[0], 'j': pair[1], 'r0_nm': r0_nm, 'label': label})

    return atom_names, atom_numbers, charges, bonds, bonds_meta

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
        if fname.endswith('xyz.xvg'):
            coord = ipt_xvg(fname)
        elif fname.endswith('f0.xvg'):
            f0 = ipt_xvg(fname)
        elif fname.endswith('f1.xvg'):
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

# def gaussian_fit(data, hist_bins=HIST_BINS):
#     """Gaussian fit to histogrammed data."""
#     N, edges = np.histogram(data, bins=hist_bins)
#     edges = edges[1:] - (edges[1] - edges[0]) / 2
#     def fun(r): return r[0]*np.exp(-((edges - r[1]) / r[2])**2) - N
#     r0 = [np.median(N), 0, 10]; lb = [0, -np.inf, 0]
#     res = least_squares(fun, r0, bounds=(lb, [np.inf, np.inf, np.inf]),
#                         xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=1_000_000)
#     fitval = res.x
#     fcurve = fitval[0]*np.exp(-((edges - fitval[1]) / fitval[2])**2)
#     return edges, N, fcurve, fitval

def gaussian_fit(data, hist_bins=HIST_BINS):
    """Gaussian fit to histogrammed data."""
    N, edges = np.histogram(data, bins=hist_bins)
    edges = edges[1:] - (edges[1] - edges[0]) / 2
    def fun(r): return r[0]*np.exp(-((edges - r[1]) / r[2])**2) + r[3]*np.exp(-((edges - r[4]) / r[5])**2) - N
    r0 = [np.max(N), -40, 10, np.max(N)/2, -20, 10]; lb = [0, -100, 1, 0, -100, 1]
    res = least_squares(fun, r0, bounds=(lb, [2*np.max(N), 0, np.inf, 2*np.max(N), 0, np.inf]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=1_000_000)
    fitval = res.x
    fcurve = fitval[0]*np.exp(-((edges - fitval[1]) / fitval[2])**2) + fitval[3]*np.exp(-((edges - fitval[4]) / fitval[5])**2)
    return edges, N, fcurve, fitval


# ================================
# Pair detection
# ================================
def pairs_from_itp_meta(bonds_meta,
                        nitrile_nm=(NITRILE_MIN_NM, NITRILE_MAX_NM),
                        carbonyl_nm=(CARBONYL_MIN_NM, CARBONYL_MAX_NM)):
    """Select nitrile/carbonyl pairs using r0 (nm) + comment label from [ bonds ]."""
    nitriles = []
    carbonyls = []
    for b in bonds_meta:
        i, j, r0, label = b['i'], b['j'], b['r0_nm'], (b['label'] or '')
        if r0 is None or not label:
            continue
        pair = (i, j) if i < j else (j, i)
        # nitrile: label says C-N and r0 within nitrile window
        if label == 'C-N' and (nitrile_nm[0] <= r0 <= nitrile_nm[1]):
            nitriles.append(pair)
        # carbonyl: label says C-O and r0 within carbonyl window
        if label == 'C-O' and (carbonyl_nm[0] <= r0 <= carbonyl_nm[1]):
            carbonyls.append(pair)
    return {'nitrile': dedup_pairs(nitriles), 'carbonyl': dedup_pairs(carbonyls)}

def guess_pairs_geometry(atom_names, atom_numbers, bonds, coord,
                         nitrile_range_A=(NITRILE_MIN_A, NITRILE_MAX_A),
                         carbonyl_range_A=(CARBONYL_MIN_A, CARBONYL_MAX_A)):
    """Fallback: guess pairs from geometry (Å) and elements."""
    avg_pos = avg_positions_from_coord(coord)
    bond_len_A = distances_for_bonds(avg_pos, bonds)
    deg = {n: 0 for n in atom_numbers}
    for a, b in bonds:
        deg[a] += 1; deg[b] += 1
    elem = {atom_numbers[i]: element_from_name(atom_names[i]) for i in range(len(atom_numbers))}

    nitriles, carbonyls = [], []
    for (a, b), dA in bond_len_A.items():
        ea, eb = elem.get(a, '?'), elem.get(b, '?')
        pair = (a, b) if a < b else (b, a)
        # nitrile candidate: C-N distance in window, and one atom terminal
        if {ea, eb} == {'C', 'N'} and nitrile_range_A[0] <= dA <= nitrile_range_A[1]:
            if deg[a] == 1 or deg[b] == 1:
                nitriles.append(pair)
        # carbonyl candidate: C-O distance in window
        if {ea, eb} == {'C', 'O'} and carbonyl_range_A[0] <= dA <= carbonyl_range_A[1]:
            carbonyls.append(pair)

    # Fallbacks to the closest if still empty
    if not nitriles:
        cn = [((a,b), dA) for (a,b), dA in bond_len_A.items()
              if {elem.get(a,'?'), elem.get(b,'?')} == {'C','N'}]
        if cn:
            nitriles.append(min(cn, key=lambda x: x[1])[0])
    if not carbonyls:
        co = [((a,b), dA) for (a,b), dA in bond_len_A.items()
              if {elem.get(a,'?'), elem.get(b,'?')} == {'C','O'}]
        if co:
            carbonyls.append(min(co, key=lambda x: x[1])[0])

    return {'nitrile': dedup_pairs(nitriles), 'carbonyl': dedup_pairs(carbonyls)}

# ================================
# Main
# ================================
def main():
    parser = argparse.ArgumentParser(description="Electric field analysis with robust nitrile/carbonyl detection from .itp bonds.")
    parser.add_argument("--pairs", nargs="*", default=None,
                        help="Atom pairs to analyze, e.g., 7,8 10,15. If omitted, auto-detect pairs.")
    parser.add_argument("--hist-bins", type=int, default=HIST_BINS, help="Histogram bins (default 128).")
    parser.add_argument("--save-prefix", default="electric_field", help="Prefix for output files.")
    # Allow customizing detection windows (Å for geometry, nm for .itp)
    parser.add_argument("--nitrile-range-A", default=f"{NITRILE_MIN_A},{NITRILE_MAX_A}",
                        help="Å range for nitrile (geometry fallback), e.g., 1.10,1.22")
    parser.add_argument("--carbonyl-range-A", default=f"{CARBONYL_MIN_A},{CARBONYL_MAX_A}",
                        help="Å range for carbonyl (geometry fallback), e.g., 1.18,1.32")
    parser.add_argument("--nitrile-range-nm", default=f"{NITRILE_MIN_NM:.6f},{NITRILE_MAX_NM:.6f}",
                        help="nm range for nitrile (from .itp bonds), e.g., 0.110000,0.122000")
    parser.add_argument("--carbonyl-range-nm", default=f"{CARBONYL_MIN_NM:.6f},{CARBONYL_MAX_NM:.6f}",
                        help="nm range for carbonyl (from .itp bonds), e.g., 0.118000,0.132000")
    args = parser.parse_args()

    nitrile_range_A = tuple(map(float, args.nitrile_range_A.split(',')))
    carbonyl_range_A = tuple(map(float, args.carbonyl_range_A.split(',')))
    nitrile_range_nm = tuple(map(float, args.nitrile_range_nm.split(',')))
    carbonyl_range_nm = tuple(map(float, args.carbonyl_range_nm.split(',')))

    # Load inputs
    itp_files = [f for f in os.listdir() if f.endswith('.itp')]
    if not itp_files:
        raise FileNotFoundError("No .itp file found in current directory.")
    atom_names, atom_idx, chg, bonds, bonds_meta = read_itp(itp_files[0])
    coord, f0, f1 = load_xvg_files()

    # Decide pairs
    if args.pairs:
        pairs = [tuple(map(int, p.split(","))) for p in args.pairs]
        source = "manual (--pairs)"
    else:
        # 1) Prefer .itp-based detection (label + r0 in nm)
        picked = pairs_from_itp_meta(bonds_meta,
                                     nitrile_nm=nitrile_range_nm,
                                     carbonyl_nm=carbonyl_range_nm)
        pairs = picked['nitrile'] + picked['carbonyl']
        source = "itp-label"
        # 2) Fallback: geometry heuristics
        if not pairs:
            picked = guess_pairs_geometry(atom_names, atom_idx, bonds, coord,
                                          nitrile_range_A=nitrile_range_A,
                                          carbonyl_range_A=carbonyl_range_A)
            pairs = picked['nitrile'] + picked['carbonyl']
            source = "geometry-fallback"
        if not pairs:
            raise RuntimeError("Auto-detection found no nitrile or carbonyl pairs. "
                               "Try adjusting detection windows or pass --pairs manually.")

    # Simple report of what we found
    print(f"[Auto-detect] Selected {len(pairs)} pair(s) via {source}: {pairs}")
    print(f"Processed {len(pairs)} pairs.")

    # Prepare plotting
    plt.figure(figsize=(5 * len(pairs), 4))
    csv_rows = []
    
    # # Get current solvent
    # current_dir = os.getcwd()
    # current_solvent = current_dir.split('/')[-2] # [-1] is 06_final_results
    
    for i, pair in enumerate(pairs, 1):
        F = compute_fields(coord, f0, f1, atom_idx, chg, pair)
        edges, N, fcurve, fitval = gaussian_fit(F, hist_bins=args.hist_bins)

        # Safe atom labels (if indices exceed table length, print indices)
        def label(atom_no):
            return atom_names[atom_no-1] if 1 <= atom_no <= len(atom_names) else f"Atom{atom_no}"
        label_pair = f"{label(pair[0])}-{label(pair[1])} ({pair[0]}-{pair[1]})"

        plt.subplot(1, len(pairs), i)
        plt.plot(edges, N, 'r.')
        plt.plot(edges, fcurve, 'b')
        meanfield = (fitval[0]*fitval[1]+fitval[3]*fitval[4])/(fitval[0]+fitval[3])
        plt.title(f"{label_pair}\n{meanfield:.2f} ± {fitval[2]/np.sqrt(2):.2f} MV/cm")
        plt.xlabel('Electric field (MV/cm)')
        plt.ylabel('Counts')

        csv_rows.append([
            label_pair, pair[0], pair[1], meanfield,
            fitval[0], fitval[1], fitval[2], fitval[2]/np.sqrt(2),
            fitval[3], fitval[4], fitval[5], fitval[5]/np.sqrt(2)
        ])

    plt.tight_layout()
    fig_path = f"{args.save_prefix}_histograms.png"
    plt.savefig(fig_path, dpi=300)
    plt.show()
    print(f"Saved plot: {fig_path}")

    # # Save individual CSV
    # csv_path = f"{args.save_prefix}_fit_results.csv"
    # with open(csv_path, "w", newline="") as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow([
    #         str(current_solvent),"Pair(Label)", "Atom1", "Atom2",
    #         "Amplitude", "Center (MV/cm)", "Sigma (MV/cm)", "Sigma/sqrt(2) (MV/cm)"
    #     ])
    #     writer.writerows(csv_rows)
    #print(f"Saved CSV:  {csv_path}")
    
    return csv_rows

if __name__ == "__main__":
    # Get subfolders of current folder
    working_dir = os.getcwd()
    print("Working directory:", str(working_dir))
    
    # batch_csv_rows = []
    # batch_csv_rows.append(["Solvent", "Pair(Label)", "Atom1", "Atom2",
    # "Amplitude", "Center (MV/cm)", "Sigma (MV/cm)", "Sigma/sqrt(2) (MV/cm)"])
    
    # Try pandas instead
    rows = []
    columns = ["Solvent", "Pair(Label)", "Atom1", "Atom2", "Meanfield",
    "Amplitude1", "Center1 (MV/cm)", "Sigma1 (MV/cm)", "Sigma1/sqrt(2) (MV/cm)",
    "Amplitude2", "Center2 (MV/cm)", "Sigma2 (MV/cm)", "Sigma2/sqrt(2) (MV/cm)"]
    
    #for folder, dirs, files in os.walk(working_dir):
    for folder in os.listdir(working_dir): # loop through folders
    
        # Look for folders with results in them    
        guessed_results = str(folder)+"/06_final_results"
        if os.path.isdir(guessed_results):
            
            # Go to folder
            os.chdir(guessed_results)

            #print(guessed_results)
            csv_rows = main()
            #batch_csv_rows.append(csv_rows)
            
            for row in csv_rows:
                clean_row = [folder] + [x.item() if isinstance(x, (np.generic,)) else x for x in row]
                rows.append(clean_row)
            
            
            # Return to working directory
            os.chdir(working_dir)
    
    # Save batch CSV after looping through folders
    csv_path = "Electric_field_batch_results.csv"
    
    # # Save CSV
    # with open(csv_path, "w", newline="") as csvfile:
    #     writer = csv.writer(csvfile)
    #     # writer.writerow([
    #     #     "Solvent", "Pair(Label)", "Atom1", "Atom2",
    #     #     "Amplitude", "Center (MV/cm)", "Sigma (MV/cm)", "Sigma/sqrt(2) (MV/cm)"
    #     # ])
    #     writer.writerows(batch_csv_rows)
        
    # Save with pandas instead
    df = pd.DataFrame(rows, columns = columns)
    df.to_csv(csv_path, index=False)
