#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 17:13:26 2025

@author: pkocheril
"""

# run this after editing the solvent topology files

targetbox = 4 # nm, edge length of desired box

from pathlib import Path
import pandas as pd
import re
import subprocess

GMX_PATH = "/groups/WeiLab/software/gromacs/bin/gmx"

def run_cmd(cmd, cwd=None):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=cwd)

def extract_section(lines, target_sections):
    """
    Extract specified sections from a list of lines.
    Returns a dict of section name -> lines and the remaining lines.
    """
    section_lines = {}
    remaining_lines = []
    current_section = None
    buffer = []

    def store_section(name, buf):
        if name in target_sections:
            section_lines[name] = list(buf)
        else:
            remaining_lines.extend(buf)

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('[') and stripped.endswith(']'):
            if current_section:
                store_section(current_section, buffer)
            current_section = stripped
            buffer = [line]
        else:
            buffer.append(line)

    if current_section:
        store_section(current_section, buffer)

    return section_lines, remaining_lines

def remove_duplicate_atomtypes(filename, output_filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    inside_atomtypes = False
    seen_atomtypes = set()
    output_lines = []

    for line in lines:
        stripped = line.strip()

        # Check for start and end of the [ atomtypes ] section
        if stripped.startswith('[ atomtypes ]'):
            inside_atomtypes = True
            output_lines.append(line)
            continue
        elif stripped.startswith('#') and inside_atomtypes:
            inside_atomtypes = False  # exiting section

        if inside_atomtypes:
            if stripped.startswith(';') or stripped == '':
                output_lines.append(line)
                continue

            atomtype = stripped.split()[0]  # first column: atom type name

            if atomtype not in seen_atomtypes:
                seen_atomtypes.add(atomtype)
                output_lines.append(line)
            else:
                # Duplicate atomtype, skip it
                continue
        else:
            output_lines.append(line)

    with open(output_filename, 'w') as f:
        f.writelines(output_lines)

    print(f"Cleaned file written to {output_filename}")

# Directory containing the .top files and solvent property table
directory = Path(".")
csv_path = directory / "Solvent_properties.csv"

# Find all XXX_###.top files and extract (XXX, ###) info
top_files = list(directory.glob("*.top"))
gro_files = list(directory.glob("*.gro"))
solvent_info = []

for file in top_files:
    match = re.match(r"([A-Za-z]{3})_(\d+)\.top", file.name)
    if match:
        solvent = match.group(1)
        volume = int(match.group(2))
        solvent_info.append({"file": file, "solvent": solvent, "volume": volume})

# Convert to DataFrame
df = pd.DataFrame(solvent_info)

# Scale volumes so they sum to 100
total_volume = df["volume"].sum()
if total_volume != 100:
    scaling_factor = 100 / total_volume
    df["scaled_volume"] = df["volume"] * scaling_factor
else:
    df["scaled_volume"] = df["volume"]

# Load solvent property table
solvent_props = pd.read_csv(csv_path)

# Check for required columns (column 2, 4, and 7 are used — index 1, 3, 6)
required_cols = [1, 3, 6]
if not all(col < solvent_props.shape[1] for col in required_cols):
    raise ValueError("CSV file does not contain the required columns (2, 4, and 7)")

# Rename relevant columns for clarity
solvent_props.columns = [f"col{i+1}" for i in range(solvent_props.shape[1])]
solvent_props = solvent_props.rename(columns={"col1": "fullname", "col2": "solvent", "col3": "formula", 
                                              "col4": "molar_mass", "col5": "BP", "col6": "MP", 
                                              "col7": "density", "col8": "solubility", "col9": "dielectric", 
                                              "col10": "flash_point", "col11": "molarity"})

# Merge dataframes on solvent abbreviation
merged = df.merge(solvent_props[["solvent", "molar_mass", "density"]], on="solvent", how="left")

# Calculate mole ratios: scaled_volume * molar_mass / density
merged["mole_ratio"] = merged["scaled_volume"] * merged["density"] / merged["molar_mass"]

merged[["file", "solvent", "volume", "scaled_volume", "density", "molar_mass", "mole_ratio"]]

targetvol = targetbox**3 # nm^3
targetvol = 1e-21*targetvol # mL
merged["targetmol"] = targetvol * 6.0221408e23 * (merged["scaled_volume"] / 100) * merged["density"] / merged["molar_mass"]

for i in range(0, merged.shape[0]):
    rawvol = merged["volume"].loc[i]
    if rawvol == 100:
        volstring = str(merged["volume"].loc[i])
    elif rawvol > 9.5:
        volstring = "0" + str(merged["volume"].loc[i])
    else:
        volstring = "00" + str(merged["volume"].loc[i])
    
    if i + 1 == merged.shape[0]: # final iteration
        boxname = "solv.gro"
    else:
        boxname = f"box{i}.gro"
    
    if i == 0:
        run_cmd([GMX_PATH, "insert-molecules", "-ci", str(merged["solvent"].loc[i] + "_" + volstring + ".gro"), "-nmol", 
                 str(round(merged["targetmol"].loc[i])), "-rot", "xyz", "-box", 
                 str(targetbox), str(targetbox), str(targetbox), "-o", boxname, 
                 "-scale", "0.4"])
    else:
        run_cmd([GMX_PATH, "insert-molecules", "-f",f"box{i-1}.gro","-ci", str(merged["solvent"].loc[i] + "_" + volstring + ".gro"), "-nmol", 
                 str(round(merged["targetmol"].loc[i])), "-rot", "xyz", "-box", 
                 str(targetbox), str(targetbox), str(targetbox), "-o", boxname, 
                 "-scale", "0.4"])
    
    # Combine .top files
    with open(str(merged["solvent"].loc[i] + "_" + volstring + ".top"), 'r') as f:
        top_lines = f.readlines()
    
    defaults_section, _ = extract_section(top_lines, ['[ defaults ]'])
    atomtypes_section, _ = extract_section(top_lines, ['[ atomtypes ]'])
    system_section, _ = extract_section(top_lines, ['[ system ]'])
    molecules_section, _ = extract_section(top_lines, ['[ molecules ]'])
    
    # Modify sections
    system_section.get('[ system ]', [])[-1] = "Solvent box\n\n"
    molecules_section.get('[ molecules ]', [])[-1] = str(merged["solvent"].loc[i] + "         " + str(round(merged["targetmol"].loc[i])))

    print(f"{i}")
    
    if i == 0:
        with open("solv.top", 'w') as f:
            f.writelines(defaults_section.get('[ defaults ]', []))
            f.writelines(atomtypes_section.get('[ atomtypes ]', []))
            f.writelines(system_section.get('[ system ]', []))
            f.writelines(molecules_section.get('[ molecules ]', []))
    else:
        with open("solv.top", 'r') as g:
            solv_lines = g.readlines()
        
        solvdefaults_section, _ = extract_section(solv_lines, ['[ defaults ]'])
        solvatomtypes_section, _ = extract_section(solv_lines, ['[ atomtypes ]'])
        solvsystem_section, _ = extract_section(solv_lines, ['[ system ]'])
        solvmolecules_section, _ = extract_section(solv_lines, ['[ molecules ]'])
        
        with open("solv.top", 'w') as f:
            f.writelines(solvdefaults_section.get('[ defaults ]', [])) # original defaults
            f.writelines(solvatomtypes_section.get('[ atomtypes ]', [])[:-1]) # original atomtypes
            # Append new atomtypes
            f.writelines(atomtypes_section.get('[ atomtypes ]', [])[2:])
            f.writelines(solvsystem_section.get('[ system ]', [])) # original system
            f.writelines(solvmolecules_section.get('[ molecules ]', [])) # original molecules
            # Append new molecules
            f.write("\n")
            f.writelines(molecules_section.get('[ molecules ]', [])[2:])

# Add .itp lines
with open("solv.top", 'r') as g:
    solv_lines = g.readlines()
solvdefaults_section, _ = extract_section(solv_lines, ['[ defaults ]'])
solvatomtypes_section, _ = extract_section(solv_lines, ['[ atomtypes ]'])
solvsystem_section, _ = extract_section(solv_lines, ['[ system ]'])
solvmolecules_section, _ = extract_section(solv_lines, ['[ molecules ]'])

with open("solv.top", 'w') as f:
    f.writelines(solvdefaults_section.get('[ defaults ]', [])) # original
    f.writelines(solvatomtypes_section.get('[ atomtypes ]', [])) # original
    #f.writelines("\n")
    # Insert #include.itp lines
    for i in range(0, merged.shape[0]):
        rawvol = merged["volume"].loc[i]
        if rawvol == 100:
            volstring = str(merged["volume"].loc[i])
        elif rawvol > 9.5:
            volstring = "0" + str(merged["volume"].loc[i])
        else:
            volstring = "00" + str(merged["volume"].loc[i])
        f.write("#include " + '"' + str(merged["solvent"].loc[i] + "_" + volstring + ".itp") + '"' + "\n")
        
    f.writelines("\n")   
    f.writelines(solvsystem_section.get('[ system ]', [])) # original system
    f.writelines(solvmolecules_section.get('[ molecules ]', [])) # original molecules

# Check for duplicate atomtypes
remove_duplicate_atomtypes("solv.top", "solv.top")
