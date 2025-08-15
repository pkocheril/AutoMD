#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:53:42 2025

@author: pkocheril
"""

top_file = "Molecule.top"
itp_file = "Molecule.itp"

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


def main(top_file, itp_file):
    # === Step 1: Extract [ defaults ] from .top ===
    with open(top_file, 'r') as f:
        top_lines = f.readlines()

    defaults_section, _ = extract_section(top_lines, ['[ defaults ]'])
    trimmed_defaults = defaults_section.get('[ defaults ]', [])[:-2]  # remove last 2 lines

    # === Step 2: Extract and CUT [ atomtypes ] from .itp ===
    with open(itp_file, 'r') as f:
        itp_lines = f.readlines()

    atomtypes_section, remaining_itp = extract_section(itp_lines, ['[ atomtypes ]'])

    with open(itp_file, 'w') as f:
        f.writelines(remaining_itp)  # Overwrite without [ atomtypes ]

    # === Step 3: Write sys.itp ===
    with open("sys.itp", 'w') as f:
        f.write("; Extracted [ defaults ] and [ atomtypes ]\n\n")
        f.writelines(trimmed_defaults)
        f.write("\n")
        f.writelines(atomtypes_section.get('[ atomtypes ]', []))

    # === Step 4: Extract [ system ] and [ molecules ] from .top ===
    system_molecules_sections, _ = extract_section(top_lines, ['[ system ]', '[ molecules ]'])

    # === Step 5: Write temp.top ===
    with open("temp.top", 'w') as f:
        f.write("; Extracted [ system ] and [ molecules ]\n\n")
        for sec in ['[ system ]', '[ molecules ]']:
            f.writelines(system_molecules_sections.get(sec, []))
            f.write("\n")

    print(f"✔️ sys.itp and temp.top generated; {itp_file} updated.")


if __name__ == "__main__":
    main(top_file, itp_file)
