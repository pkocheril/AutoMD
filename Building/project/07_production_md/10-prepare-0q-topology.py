#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 09:26:19 2025

@author: pkocheril
"""

topoltop = "topol.top"
topol0q = "topol_0q.top"

import re

def create_topol_0q(topol_file, output_file):
    """
    Copy topol_file to output_file, modifying #include lines to use
    the _0q.itp copies instead of the original .itp files.
    """
    include_pattern = re.compile(r'#include\s+"itp/([A-Za-z]+_\d+)\.itp"')

    with open(topol_file, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        match = include_pattern.search(line)
        if match:
            base_name = match.group(1)
            new_line = f'#include "itp/{base_name}_0q.itp"\n'
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(new_lines)

    print(f"Created {output_file} with updated #include lines")

# Make 0q topology
create_topol_0q(topoltop, topol0q)