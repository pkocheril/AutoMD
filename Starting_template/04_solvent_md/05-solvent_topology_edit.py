#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 15:40:31 2025

@author: pkocheril
"""

import os
import re
import shutil

def remove_unwanted_sections(lines):
    """
    Remove specific sections ([defaults], [atomtypes], [system], [molecules]) from a GROMACS file.
    """
    unwanted = {'defaults', 'atomtypes', 'system', 'molecules'}
    cleaned_lines = []
    keep = True

    for line in lines:
        section_match = re.match(r'\[\s*(\w+)\s*\]', line.strip().lower())
        if section_match:
            section = section_match.group(1)
            keep = section not in unwanted

        if keep:
            cleaned_lines.append(line)

    return cleaned_lines


def fix_charge(lines):
    inside_atoms = False
    fixed_lines = []
    atom_lines = []
    parsed_atom_lines = []
    original_atom_lines = []
    after_atoms = []

    for i, line in enumerate(lines):
        if line.strip().lower().startswith('[ atoms ]'):
            inside_atoms = True
            fixed_lines.append(line)
            continue
        if inside_atoms:
            if re.match(r'\[\s*\w+\s*\]', line) or i == len(lines) - 1:
                inside_atoms = False
                if i == len(lines) - 1:  # include last line
                    atom_lines.append(line)
                after_atoms = lines[i:]  # preserve rest
                break
            atom_lines.append(line)

    # Parse atom lines
    for line in atom_lines:
        original_atom_lines.append(line)  # Save for later output
        if line.strip() == '' or line.strip().startswith(';'):
            parsed_atom_lines.append((line, None))
            continue
        parts = line.split()
        if len(parts) < 7:
            parsed_atom_lines.append((line, None))
            continue
        try:
            charge = float(parts[6])
            parsed_atom_lines.append((parts.copy(), charge))  # copy to avoid overwriting
        except ValueError:
            parsed_atom_lines.append((line, None))

    # Sum charges
    charges = [x[1] for x in parsed_atom_lines if x[1] is not None]
    total_charge = sum(charges)

    # Fix charge if needed
    if abs(total_charge) > 1e-6:
        for i in reversed(range(len(parsed_atom_lines))):
            item, charge = parsed_atom_lines[i]
            if charge is not None:
                corrected = charge - total_charge
                item[6] = f"{corrected:.6f}"
                parsed_atom_lines[i] = (item, corrected)
                print(f"⚠️ Corrected charge by {(-total_charge):+.6f} in atom line {i+1}")
                break

    # Output original lines, commented
    fixed_lines.append("; Original atom lines:\n")
    for line in original_atom_lines:
        if line.strip() == '':
            fixed_lines.append(";\n")
        elif line.strip().startswith(';'):
            fixed_lines.append(line)
        else:
            fixed_lines.append("; " + line if not line.startswith(';') else line)

    # Output corrected atom lines
    fixed_lines.append("; Corrected charges below\n")
    for item, charge in parsed_atom_lines:
        if isinstance(item, list):
            formatted = "{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:>10} {:>10}".format(
                *item[:7], item[7] if len(item) > 7 else ''
            )
            fixed_lines.append(formatted + '\n')

    # Append remainder of file
    fixed_lines.append('\n')
    fixed_lines.extend(after_atoms)
    return fixed_lines

def replace_mol(lines, new_name):
    return [line.replace("MOL", new_name) for line in lines]

def process_file(top_file):
    base_name = os.path.splitext(top_file)[0]
    match = re.match(r'([A-Za-z]{3})_\d{1,3}', base_name)
    if not match:
        print(f"⚠️ Skipping {top_file}: filename does not match expected pattern.")
        return
    solvent = match.group(1)

    itp_file = f"{base_name}.itp"
    shutil.copyfile(top_file, itp_file)

    with open(top_file, 'r') as f:
        lines = f.readlines()

    lines = remove_unwanted_sections(lines)

    # Fix net charge
    lines = fix_charge(lines)

    # Add [ moleculetype ] section manually at the top
    moleculetype_block = [
        "[ moleculetype ]\n",
        ";name            nrexcl\n",
        " MOL              3\n\n"
    ]
    
    # Prepend it to the cleaned lines
    lines = moleculetype_block + lines

    # Replace "MOL" with solvent name throughout
    lines = replace_mol(lines, solvent)

    with open(itp_file, 'w') as f:
        f.writelines(lines)

    # Update matching .gro file
    gro_file = f"{base_name}.gro"
    if os.path.exists(gro_file):
        with open(gro_file, 'r') as f:
            gro_lines = f.readlines()
        gro_lines = replace_mol(gro_lines, solvent)
        with open(gro_file, 'w') as f:
            f.writelines(gro_lines)

    print(f"✅ Processed: {top_file} → {itp_file} + charge check + label rename")

if __name__ == "__main__":
    for file in os.listdir():
        if file.endswith(".top"):
            process_file(file)

