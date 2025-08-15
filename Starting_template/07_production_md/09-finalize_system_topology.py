#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 15:16:37 2025

@author: pkocheril
"""

itpfolder = "itp"
solventtop = f"{itpfolder}/solv.top"
sysitp = f"{itpfolder}/sys.itp"
sysmerge = f"{itpfolder}/sys_merged.itp"
sysclean = f"{itpfolder}/sys_clean.itp"
temptop = "temp.top"
topoltop = "topol.top"
edgedist = "2.0" # nm, half of box edge length
topol0q = "topol_0q.top"

import os
import re
import shutil
# import subprocess

# GMX_PATH = "/groups/WeiLab/software/gromacs/bin/gmx"

# def run_cmd(cmd, cwd=None):
#     print(f"Running: {' '.join(cmd)}")
#     subprocess.run(cmd, check=True, cwd=cwd)

def append_atomtypes(source_file, target_file, output_file):
    """
    Append [ atomtypes ] entries from source_file into target_file
    directly after the header line in target_file, removing blank lines,
    skipping '#include' lines, and avoiding duplicate section headers.
    """
    # Read files
    with open(source_file, 'r') as f:
        source_lines = f.readlines()
    with open(target_file, 'r') as f:
        target_lines = f.readlines()

    # --- Extract entries from source ---
    in_section = False
    atomtypes_to_add = []
    for line in source_lines:
        stripped = line.strip()
        if stripped.startswith('[ atomtypes ]'):
            in_section = True
            continue  # skip header
        elif stripped.startswith('[') and ']' in stripped and in_section:
            break  # end of section
        elif in_section:
            # Skip '#include' lines
            if stripped.lower().startswith('#include'):
                continue
            atomtypes_to_add.append(line)

    # Remove leading blanks/comments from source section
    while atomtypes_to_add and (atomtypes_to_add[0].strip() == '' or atomtypes_to_add[0].strip().startswith(';')):
        atomtypes_to_add.pop(0)

    # --- Merge into target ---
    merged_lines = []
    in_section = False
    for i, line in enumerate(target_lines):
        merged_lines.append(line)
        if line.strip().startswith('[ atomtypes ]') and not in_section:
            in_section = True
            # Find first blank line after the header
            j = i + 1
            while j < len(target_lines) and target_lines[j].strip() != '':
                merged_lines.append(target_lines[j])
                j += 1
            # Append new entries right after the original block
            merged_lines.extend(atomtypes_to_add)
            merged_lines.append("\n")  # keep one blank line after appended entries
            # Skip over the original blank line(s)
            while j < len(target_lines) and target_lines[j].strip() == '':
                j += 1
            # Resume copying rest of file
            merged_lines.extend(target_lines[j:])
            break

    # --- Write merged file ---
    with open(output_file, 'w') as f:
        f.writelines(merged_lines)

    print(f"[ atomtypes ] from {source_file} appended to {target_file} → {output_file}")

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

def insert_itp_includes_before_system(topol_file, itp_folder, output_file):
    """
    Insert #include lines for all .itp files in itp_folder into topol_file,
    excluding sys.itp and sys_merged.itp, placing them right before the [ system ] section.
    Avoids adding duplicates if include lines are already present.
    Files containing 'sys' appear first in the list.
    """
    # Gather .itp files
    itp_files = [
        f for f in os.listdir(itp_folder)
        if f.endswith(".itp") and f not in ("sys.itp", "sys_merged.itp")
    ]
    # Sort so 'sys' files come first, each group sorted alphabetically
    itp_files.sort(key=lambda x: (0 if "sys" in x.lower() else 1, x.lower()))

    # Read topol.top
    with open(topol_file, 'r') as f:
        topol_lines = f.readlines()

    # Build list of include lines to add
    existing_lines = {line.strip() for line in topol_lines}
    include_lines = [
        f'#include "itp/{fname}"\n'
        for fname in itp_files
        if f'#include "itp/{fname}"' not in existing_lines
    ]

    if not include_lines:
        print("No new includes to add — all already present.")
        with open(output_file, 'w') as f:
            f.writelines(topol_lines)
        return

    # Insert before [ system ]
    merged_lines = []
    inserted = False
    for line in topol_lines:
        stripped = line.strip()
        if stripped.lower().startswith("[ system ]") and not inserted:
            if merged_lines and not merged_lines[-1].endswith("\n"):
                merged_lines[-1] += "\n"
            merged_lines.extend(include_lines)
            merged_lines.append("\n")
            inserted = True
        merged_lines.append(line)

    # Write updated file
    with open(output_file, 'w') as f:
        f.writelines(merged_lines)

    print(f"Inserted {len(include_lines)} new include lines before [ system ] in {output_file}")

# def insert_itp_includes_before_system(topol_file, itp_folder, output_file):
#     """
#     Insert #include lines for all .itp files in itp_folder into topol_file,
#     excluding sys.itp and sys_merged.itp, placing them right before the [ system ] section.
#     Avoids adding duplicates if include lines are already present.
#     """
#     # Gather .itp files
#     itp_files = [
#         f for f in os.listdir(itp_folder)
#         if f.endswith(".itp") and f not in ("sys.itp", "sys_merged.itp")
#     ]
#     itp_files.sort()  # alphabetical order for consistency

#     # Read topol.top
#     with open(topol_file, 'r') as f:
#         topol_lines = f.readlines()

#     # Build list of include lines to add
#     existing_lines = {line.strip() for line in topol_lines}
#     include_lines = [
#         f'#include "itp/{fname}"\n'
#         for fname in itp_files
#         if f'#include "itp/{fname}"' not in existing_lines
#     ]

#     if not include_lines:
#         print("No new includes to add — all already present.")
#         with open(output_file, 'w') as f:
#             f.writelines(topol_lines)
#         return

#     # Insert before [ system ]
#     merged_lines = []
#     inserted = False
#     for line in topol_lines:
#         stripped = line.strip()
#         if stripped.lower().startswith("[ system ]") and not inserted:
#             # Ensure a newline before the inserted lines
#             if merged_lines and not merged_lines[-1].endswith("\n"):
#                 merged_lines[-1] += "\n"
#             merged_lines.extend(include_lines)
#             merged_lines.append("\n")  # blank line before [ system ]
#             inserted = True
#         merged_lines.append(line)

#     # Write updated file
#     with open(output_file, 'w') as f:
#         f.writelines(merged_lines)

#     print(f"Inserted {len(include_lines)} new include lines before [ system ] in {output_file}")

def make_zero_charge_copies(itp_folder):
    """
    For each file matching XXX_###.itp in itp_folder, create a copy named XXX_###_0q.itp.
    In the copy, set all charges in the [ atoms ] section (column 7) to 0.000000.
    """
    pattern = re.compile(r"^[A-Za-z]+_\d+\.itp$")
    files_to_copy = [f for f in os.listdir(itp_folder) if pattern.match(f)]

    for filename in files_to_copy:
        src_path = os.path.join(itp_folder, filename)
        new_filename = filename.replace(".itp", "_0q.itp")
        dst_path = os.path.join(itp_folder, new_filename)

        # Copy the file
        shutil.copy(src_path, dst_path)

        # Edit charges in the new file
        with open(dst_path, 'r') as f:
            lines = f.readlines()

        output_lines = []
        inside_atoms = False
        for line in lines:
            stripped = line.strip()

            if stripped.startswith("[ atoms ]"):
                inside_atoms = True
                output_lines.append(line)
                continue
            elif stripped.startswith("[") and "]" in stripped and inside_atoms:
                inside_atoms = False  # leaving section

            if inside_atoms:
                # Skip comments or blank lines
                if stripped.startswith(";") or stripped == "":
                    output_lines.append(line)
                    continue

                cols = stripped.split()
                if len(cols) >= 7:
                    cols[6] = "0.000000"
                    line = " ".join(cols) + "\n"
            output_lines.append(line)

        with open(dst_path, 'w') as f:
            f.writelines(output_lines)

        print(f"Created zero-charge copy: {new_filename}")



########### Run functions

# Copy atomtypes to sys.itp
append_atomtypes(solventtop, sysitp, sysmerge)
remove_duplicate_atomtypes(sysmerge, sysclean)

# Put #include itp lines before [ system ] in topol.top
insert_itp_includes_before_system(temptop, itpfolder, topoltop)

# Make 0q files
make_zero_charge_copies(itpfolder)

# # Put molecule inside a box
# molgro = [
#     f for f in os.listdir()
#     if f.endswith(".gro") and f not in ("solv_p.gro")
# ]
# run_cmd([GMX_PATH, "editconf", "-f", molgro, "-o", "box0.gro", "-c", "-d", edgedist, "-bt", "cubic"])

# # Solvate
# run_cmd([GMX_PATH, "solvate", "-cp", "box0.gro", "-cs", "solv_p.gro", "-o", "MIX.gro", "-p", topoltop])






# def append_atomtypes(source_file, target_file, output_file):
#     """
#     Append the [ atomtypes ] section from source_file to the [ atomtypes ] section in target_file.
#     """
#     # Read source and target files
#     with open(source_file, 'r') as f:
#         source_lines = f.readlines()
#     with open(target_file, 'r') as f:
#         target_lines = f.readlines()

#     # --- Extract [ atomtypes ] from source ---
#     in_section = False
#     atomtypes_to_add = []
#     for line in source_lines:
#         stripped = line.strip()
#         if stripped.startswith('[ atomtypes ]'):
#             in_section = True
#             atomtypes_to_add.append(line)  # keep the header
#             continue
#         elif stripped.startswith('[') and ']' in stripped and in_section:
#             break  # end of section
#         elif in_section:
#             atomtypes_to_add.append(line)

#     # --- Append to target's [ atomtypes ] ---
#     in_section = False
#     merged_lines = []
#     appended = False
#     for line in target_lines:
#         stripped = line.strip()
#         if stripped.startswith('[ atomtypes ]'):
#             in_section = True
#             merged_lines.append(line)
#             continue
#         elif stripped.startswith('[') and ']' in stripped and in_section:
#             # End of atomtypes section
#             merged_lines.extend(atomtypes_to_add[1:])  # skip header in appended part
#             in_section = False
#             appended = True
#         merged_lines.append(line)

#     # If we never hit a new section, append at the end
#     if in_section and not appended:
#         merged_lines.extend(atomtypes_to_add[1:])

#     # --- Write merged file ---
#     with open(output_file, 'w') as f:
#         f.writelines(merged_lines)

#     print(f"[ atomtypes ] from {source_file} appended to {target_file} → {output_file}")



# code from 08-build_box, maybe some useful functions there

# targetbox = 4 # nm, edge length of desired box

# from pathlib import Path
# import pandas as pd
# import re
# import subprocess

# GMX_PATH = "/groups/WeiLab/software/gromacs/bin/gmx"

# def run_cmd(cmd, cwd=None):
#     print(f"Running: {' '.join(cmd)}")
#     subprocess.run(cmd, check=True, cwd=cwd)

# def extract_section(lines, target_sections):
#     """
#     Extract specified sections from a list of lines.
#     Returns a dict of section name -> lines and the remaining lines.
#     """
#     section_lines = {}
#     remaining_lines = []
#     current_section = None
#     buffer = []

#     def store_section(name, buf):
#         if name in target_sections:
#             section_lines[name] = list(buf)
#         else:
#             remaining_lines.extend(buf)

#     for line in lines:
#         stripped = line.strip()
#         if stripped.startswith('[') and stripped.endswith(']'):
#             if current_section:
#                 store_section(current_section, buffer)
#             current_section = stripped
#             buffer = [line]
#         else:
#             buffer.append(line)

#     if current_section:
#         store_section(current_section, buffer)

#     return section_lines, remaining_lines

# def remove_duplicate_atomtypes(filename, output_filename):
#     with open(filename, 'r') as f:
#         lines = f.readlines()

#     inside_atomtypes = False
#     seen_atomtypes = set()
#     output_lines = []

#     for line in lines:
#         stripped = line.strip()

#         # Check for start and end of the [ atomtypes ] section
#         if stripped.startswith('[ atomtypes ]'):
#             inside_atomtypes = True
#             output_lines.append(line)
#             continue
#         elif stripped.startswith('#') and inside_atomtypes:
#             inside_atomtypes = False  # exiting section

#         if inside_atomtypes:
#             if stripped.startswith(';') or stripped == '':
#                 output_lines.append(line)
#                 continue

#             atomtype = stripped.split()[0]  # first column: atom type name

#             if atomtype not in seen_atomtypes:
#                 seen_atomtypes.add(atomtype)
#                 output_lines.append(line)
#             else:
#                 # Duplicate atomtype, skip it
#                 continue
#         else:
#             output_lines.append(line)

#     with open(output_filename, 'w') as f:
#         f.writelines(output_lines)

#     print(f"Cleaned file written to {output_filename}")
    
# # Directory containing the .top files and solvent property table
# directory = Path(".")
# csv_path = directory / "Solvent_properties.csv"

# # Find all XXX_###.top files and extract (XXX, ###) info
# top_files = list(directory.glob("*.top"))
# gro_files = list(directory.glob("*.gro"))
# solvent_info = []

# for file in top_files:
#     match = re.match(r"([A-Za-z]{3})_(\d+)\.top", file.name)
#     if match:
#         solvent = match.group(1)
#         volume = int(match.group(2))
#         solvent_info.append({"file": file, "solvent": solvent, "volume": volume})

# # Convert to DataFrame
# df = pd.DataFrame(solvent_info)

# # Scale volumes so they sum to 100
# total_volume = df["volume"].sum()
# if total_volume != 100:
#     scaling_factor = 100 / total_volume
#     df["scaled_volume"] = df["volume"] * scaling_factor
# else:
#     df["scaled_volume"] = df["volume"]

# # Load solvent property table
# solvent_props = pd.read_csv(csv_path)

# # Check for required columns (column 2, 4, and 7 are used — index 1, 3, 6)
# required_cols = [1, 3, 6]
# if not all(col < solvent_props.shape[1] for col in required_cols):
#     raise ValueError("CSV file does not contain the required columns (2, 4, and 7)")

# # Rename relevant columns for clarity
# solvent_props.columns = [f"col{i+1}" for i in range(solvent_props.shape[1])]
# solvent_props = solvent_props.rename(columns={"col1": "fullname", "col2": "solvent", "col3": "formula", 
#                                               "col4": "molar_mass", "col5": "BP", "col6": "MP", 
#                                               "col7": "density", "col8": "solubility", "col9": "dielectric", 
#                                               "col10": "flash_point", "col11": "molarity"})

# # Merge dataframes on solvent abbreviation
# merged = df.merge(solvent_props[["solvent", "molar_mass", "density"]], on="solvent", how="left")

# # Calculate mole ratios: scaled_volume * molar_mass / density
# merged["mole_ratio"] = merged["scaled_volume"] * merged["density"] / merged["molar_mass"]

# merged[["file", "solvent", "volume", "scaled_volume", "density", "molar_mass", "mole_ratio"]]

# targetvol = targetbox**3 # nm^3
# targetvol = 1e-21*targetvol # mL
# merged["targetmol"] = targetvol * 6.0221408e23 * (merged["scaled_volume"] / 100) * merged["density"] / merged["molar_mass"]

# for i in range(0, merged.shape[0]):
#     rawvol = merged["volume"].loc[i]
#     if rawvol == 100:
#         volstring = str(merged["volume"].loc[i])
#     elif rawvol > 9.5:
#         volstring = "0" + str(merged["volume"].loc[i])
#     else:
#         volstring = "00" + str(merged["volume"].loc[i])
    
#     if i + 1 == merged.shape[0]: # final iteration
#         boxname = "solv.gro"
#     else:
#         boxname = f"box{i}.gro"
    
#     if i == 0:
#         run_cmd([GMX_PATH, "insert-molecules", "-ci", str(merged["solvent"].loc[i] + "_" + volstring + ".gro"), "-nmol", 
#                  str(round(merged["targetmol"].loc[i])), "-rot", "xyz", "-box", 
#                  str(targetbox), str(targetbox), str(targetbox), "-o", boxname, 
#                  "-scale", "0.4"])
#     else:
#         run_cmd([GMX_PATH, "insert-molecules", "-f",f"box{i-1}.gro","-ci", str(merged["solvent"].loc[i] + "_" + volstring + ".gro"), "-nmol", 
#                  str(round(merged["targetmol"].loc[i])), "-rot", "xyz", "-box", 
#                  str(targetbox), str(targetbox), str(targetbox), "-o", boxname, 
#                  "-scale", "0.4"])
    
#     # Combine .top files
#     with open(str(merged["solvent"].loc[i] + "_" + volstring + ".top"), 'r') as f:
#         top_lines = f.readlines()
    
#     defaults_section, _ = extract_section(top_lines, ['[ defaults ]'])
#     atomtypes_section, _ = extract_section(top_lines, ['[ atomtypes ]'])
#     system_section, _ = extract_section(top_lines, ['[ system ]'])
#     molecules_section, _ = extract_section(top_lines, ['[ molecules ]'])
    
#     # Modify sections
#     system_section.get('[ system ]', [])[-1] = "Solvent box\n\n"
#     molecules_section.get('[ molecules ]', [])[-1] = str(merged["solvent"].loc[i] + "         " + str(round(merged["targetmol"].loc[i])))

#     print(f"{i}")
    
#     if i == 0:
#         with open("solv.top", 'w') as f:
#             f.writelines(defaults_section.get('[ defaults ]', []))
#             f.writelines(atomtypes_section.get('[ atomtypes ]', []))
#             f.writelines(system_section.get('[ system ]', []))
#             f.writelines(molecules_section.get('[ molecules ]', []))
#     else:
#         with open("solv.top", 'r') as g:
#             solv_lines = g.readlines()
        
#         solvdefaults_section, _ = extract_section(solv_lines, ['[ defaults ]'])
#         solvatomtypes_section, _ = extract_section(solv_lines, ['[ atomtypes ]'])
#         solvsystem_section, _ = extract_section(solv_lines, ['[ system ]'])
#         solvmolecules_section, _ = extract_section(solv_lines, ['[ molecules ]'])
        
#         with open("solv.top", 'w') as f:
#             f.writelines(solvdefaults_section.get('[ defaults ]', [])) # original defaults
#             f.writelines(solvatomtypes_section.get('[ atomtypes ]', [])[:-1]) # original atomtypes
#             # Append new atomtypes
#             f.writelines(atomtypes_section.get('[ atomtypes ]', [])[2:])
#             f.writelines(solvsystem_section.get('[ system ]', [])) # original system
#             f.writelines(solvmolecules_section.get('[ molecules ]', [])) # original molecules
#             # Append new molecules
#             f.write("\n")
#             f.writelines(molecules_section.get('[ molecules ]', [])[2:])

# # Add .itp lines
# with open("solv.top", 'r') as g:
#     solv_lines = g.readlines()
# solvdefaults_section, _ = extract_section(solv_lines, ['[ defaults ]'])
# solvatomtypes_section, _ = extract_section(solv_lines, ['[ atomtypes ]'])
# solvsystem_section, _ = extract_section(solv_lines, ['[ system ]'])
# solvmolecules_section, _ = extract_section(solv_lines, ['[ molecules ]'])

# with open("solv.top", 'w') as f:
#     f.writelines(solvdefaults_section.get('[ defaults ]', [])) # original
#     f.writelines(solvatomtypes_section.get('[ atomtypes ]', [])) # original
#     #f.writelines("\n")
#     # Insert #include.itp lines
#     for i in range(0, merged.shape[0]):
#         rawvol = merged["volume"].loc[i]
#         if rawvol == 100:
#             volstring = str(merged["volume"].loc[i])
#         elif rawvol > 9.5:
#             volstring = "0" + str(merged["volume"].loc[i])
#         else:
#             volstring = "00" + str(merged["volume"].loc[i])
#         f.write("#include " + '"' + str(merged["solvent"].loc[i] + "_" + volstring + ".itp") + '"' + "\n")
        
#     f.writelines("\n")   
#     f.writelines(solvsystem_section.get('[ system ]', [])) # original system
#     f.writelines(solvmolecules_section.get('[ molecules ]', [])) # original molecules

# # Check for duplicate atomtypes
# remove_duplicate_atomtypes("solv.top", "solv.top")
