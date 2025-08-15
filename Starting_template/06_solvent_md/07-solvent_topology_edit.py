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

# import os
# import re
# import shutil

# def remove_sections(lines, sections_to_remove):
#     result = []
#     skip = False
#     for line in lines:
#         if re.match(r'\[\s*(\w+)\s*\]', line):
#             section = re.findall(r'\[\s*(\w+)\s*\]', line)[0].lower()
#             skip = section in sections_to_remove
#         if not skip:
#             result.append(line)
#     return result

# def fix_charge(lines):
#     inside_atoms = False
#     fixed_lines = []
#     atom_lines = []
#     parsed_atom_lines = []
#     original_atom_lines = []
#     after_atoms = []

#     for i, line in enumerate(lines):
#         if line.strip().lower().startswith('[ atoms ]'):
#             inside_atoms = True
#             fixed_lines.append(line)
#             continue
#         if inside_atoms:
#             if re.match(r'\[\s*\w+\s*\]', line) or i == len(lines) - 1:
#                 inside_atoms = False
#                 if i == len(lines) - 1:  # include last line
#                     atom_lines.append(line)
#                 after_atoms = lines[i:]  # preserve rest
#                 break
#             atom_lines.append(line)

#     # Parse atom lines
#     for line in atom_lines:
#         original_atom_lines.append(line)  # Save for later output
#         if line.strip() == '' or line.strip().startswith(';'):
#             parsed_atom_lines.append((line, None))
#             continue
#         parts = line.split()
#         if len(parts) < 7:
#             parsed_atom_lines.append((line, None))
#             continue
#         try:
#             charge = float(parts[6])
#             parsed_atom_lines.append((parts.copy(), charge))  # copy to avoid overwriting
#         except ValueError:
#             parsed_atom_lines.append((line, None))

#     # Sum charges
#     charges = [x[1] for x in parsed_atom_lines if x[1] is not None]
#     total_charge = sum(charges)

#     # Fix charge if needed
#     if abs(total_charge) > 1e-6:
#         for i in reversed(range(len(parsed_atom_lines))):
#             item, charge = parsed_atom_lines[i]
#             if charge is not None:
#                 corrected = charge - total_charge
#                 item[6] = f"{corrected:.6f}"
#                 parsed_atom_lines[i] = (item, corrected)
#                 print(f"⚠️ Corrected charge by {(-total_charge):+.6f} in atom line {i+1}")
#                 break

#     # Output original lines, commented
#     fixed_lines.append("; Original atom lines:\n")
#     for line in original_atom_lines:
#         if line.strip() == '':
#             fixed_lines.append(";\n")
#         elif line.strip().startswith(';'):
#             fixed_lines.append(line)
#         else:
#             fixed_lines.append("; " + line if not line.startswith(';') else line)

#     # Output corrected atom lines
#     fixed_lines.append("; Corrected charges below\n")
#     for item, charge in parsed_atom_lines:
#         if isinstance(item, list):
#             formatted = "{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:>10} {:>10}".format(
#                 *item[:7], item[7] if len(item) > 7 else ''
#             )
#             fixed_lines.append(formatted + '\n')

#     # Append remainder of file
#     fixed_lines.append('\n')
#     fixed_lines.extend(after_atoms)
#     return fixed_lines

# def replace_mol(lines, new_name):
#     return [line.replace("MOL", new_name) for line in lines]

# def process_file(top_file):
#     base_name = os.path.splitext(top_file)[0]
#     match = re.match(r'([A-Za-z]{3})_\d{1,3}', base_name)
#     if not match:
#         print(f"⚠️ Skipping {top_file}: filename does not match expected pattern.")
#         return
#     solvent = match.group(1)

#     itp_file = f"{base_name}.itp"
#     shutil.copyfile(top_file, itp_file)

#     with open(itp_file, 'r') as f:
#         lines = f.readlines()

#     # Remove specific sections
#     sections_to_remove = ['defaults', 'atomtypes', 'system', 'molecules']
#     lines = remove_sections(lines, sections_to_remove)

#     # Fix net charge
#     lines = fix_charge(lines)

#     # Replace "MOL" with solvent name throughout
#     lines = replace_mol(lines, solvent)

#     with open(itp_file, 'w') as f:
#         f.writelines(lines)

#     # Update matching .gro file
#     gro_file = f"{base_name}.gro"
#     if os.path.exists(gro_file):
#         with open(gro_file, 'r') as f:
#             gro_lines = f.readlines()
#         gro_lines = replace_mol(gro_lines, solvent)
#         with open(gro_file, 'w') as f:
#             f.writelines(gro_lines)

#     print(f"✅ Processed: {top_file} → {itp_file} + charge check + label rename")

# if __name__ == "__main__":
#     for file in os.listdir():
#         if file.endswith(".top"):
#             process_file(file)


# import os
# import re
# import shutil

# def extract_sections(lines):
#     """Return a dict of section_name -> lines."""
#     sections = {}
#     current_section = None
#     buffer = []
#     for line in lines:
#         if line.strip().startswith('[') and line.strip().endswith(']'):
#             if current_section:
#                 sections[current_section] = list(buffer)
#             current_section = line.strip()
#             buffer = [line]
#         else:
#             buffer.append(line)
#     if current_section:
#         sections[current_section] = list(buffer)
#     return sections

# def remove_sections(lines, sections_to_remove):
#     """Return lines with specified sections removed."""
#     sections = extract_sections(lines)
#     cleaned_lines = []
#     for key, value in sections.items():
#         if key.lower() not in sections_to_remove:
#             cleaned_lines.extend(value)
#     return cleaned_lines

# def fix_charge(lines):
#     """Ensure net charge in [ atoms ] section is exactly zero by adjusting the last atom."""
#     inside_atoms = False
#     fixed_lines = []
#     atom_lines = []
#     pre_atom_lines = []

#     for line in lines:
#         if line.strip().lower().startswith('[ atoms ]'):
#             inside_atoms = True
#             fixed_lines.append(line)
#             continue
#         if inside_atoms:
#             if line.strip().startswith('['):
#                 inside_atoms = False
#                 break
#             atom_lines.append(line)
#         else:
#             fixed_lines.append(line)

#     # Extract and correct charges
#     charges = []
#     parsed_atom_lines = []
#     for line in atom_lines:
#         if line.strip() == '' or line.strip().startswith(';'):
#             parsed_atom_lines.append((line, None))  # Preserve as-is
#             continue
#         parts = line.split()
#         if len(parts) < 7:
#             parsed_atom_lines.append((line, None))
#             continue
#         try:
#             charge = float(parts[6])
#             charges.append(charge)
#             parsed_atom_lines.append((parts, charge))
#         except ValueError:
#             parsed_atom_lines.append((line, None))

#     # Compute correction
#     total_charge = sum(charges)
#     if abs(total_charge) > 1e-6:
#         # Apply correction to last atom line with valid charge
#         for i in reversed(range(len(parsed_atom_lines))):
#             item, charge = parsed_atom_lines[i]
#             if charge is not None:
#                 corrected_charge = charge - total_charge
#                 item[6] = f"{corrected_charge:.6f}"
#                 parsed_atom_lines[i] = (item, corrected_charge)
#                 print(f"⚠️ Adjusted charge by {(-total_charge):+.6f} to ensure neutrality.")
#                 break

#     # Rebuild atom section
#     for item, charge in parsed_atom_lines:
#         if isinstance(item, list):
#             formatted = "{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:>10}     {:<8} ; qtot {:.3f}\n".format(
#                 *item[:7], item[7] if len(item) > 7 else '', sum(float(i[1]) for i in parsed_atom_lines if i[1] is not None)
#             )
#             fixed_lines.append(formatted)
#         else:
#             fixed_lines.append(item)

#     # Append the remaining lines after [ atoms ] section
#     remaining = False
#     for line in lines:
#         if remaining:
#             fixed_lines.append(line)
#         if line.strip().startswith('[ atoms ]'):
#             remaining = True

#     return fixed_lines


# def replace_mol_with_tag(lines, tag):
#     """Replace 'MOL' with solvent tag in [ moleculetype ] section only."""
#     result = []
#     in_moleculetype = False
#     for line in lines:
#         if line.strip().lower().startswith('[ moleculetype ]'):
#             in_moleculetype = True
#             result.append(line)
#             continue
#         if in_moleculetype:
#             if line.strip() == '' or line.strip().startswith(';'):
#                 result.append(line)
#                 continue
#             parts = line.strip().split()
#             if len(parts) >= 1 and parts[0] == 'MOL':
#                 parts[0] = tag
#                 result.append('  '.join(parts) + '\n')
#             else:
#                 result.append(line)
#             in_moleculetype = False
#         else:
#             result.append(line)
#     return result

# def replace_mol_in_gro(gro_path, tag):
#     with open(gro_path, 'r') as f:
#         lines = f.readlines()
#     new_lines = [line.replace('MOL', tag) for line in lines]
#     with open(gro_path, 'w') as f:
#         f.writelines(new_lines)

# def main():
#     for filename in os.listdir():
#         if filename.endswith('.top'):
#             tag_match = re.match(r'([A-Za-z]{3})_\d{1,3}\.top$', filename)
#             if not tag_match:
#                 continue
#             tag = tag_match.group(1)
#             print(f"\n📦 Processing: {filename} | Solvent Tag: {tag}")

#             # Step 1: Copy to .itp
#             base = os.path.splitext(filename)[0]
#             itp_filename = f"{base}.itp"
#             shutil.copy(filename, itp_filename)

#             with open(itp_filename, 'r') as f:
#                 lines = f.readlines()

#             # Step 2: Remove unwanted sections
#             lines = remove_sections(lines, ['[ defaults ]', '[ atomtypes ]', '[ system ]', '[ molecules ]'])

#             # Step 3: Check/fix charge
#             lines = fix_charge(lines)

#             # Step 4: Rename MOL → XXX in moleculetype
#             lines = replace_mol_with_tag(lines, tag)

#             with open(itp_filename, 'w') as f:
#                 f.writelines(lines)
#             print(f"✅ Wrote {itp_filename}")

#     # Step 5: Modify .gro files
#     for filename in os.listdir():
#         if filename.endswith('.gro'):
#             tag_match = re.match(r'([A-Za-z]{3})_\d{1,3}\.gro$', filename)
#             if not tag_match:
#                 continue
#             tag = tag_match.group(1)
#             print(f"🧬 Fixing .gro file: {filename}")
#             replace_mol_in_gro(filename, tag)

#     print("\n✅ All done.")

# if __name__ == "__main__":
#     main()
