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
