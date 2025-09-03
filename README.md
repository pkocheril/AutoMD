# AutoMD
Automated molecular dynamics simulations on small molecules in solvents. Also includes downstream analysis tools for calculating and visualizing electric field distributions.

Written in Bash and Python (3.11.13).


## Highlights and main functions
* Automated 3D geometry building and coarse minimization with OpenBabel
* Generates a linked job file for Gaussian16, containing:
  * Rough geometry optimization (HF/STO-3G)
  * Opt+freq DFT (B3LYP/6-31G(d,p))
  * Single-point energy calculations in vacuum and water
* Automatic solute parametrization with Multiwfn and Sobtop (RESP2)
* Building and equilibrating of a solvent box from bulk solvent properties (density and molecular weight)
* Production MD with a solvated molecule
* Calculation of electric fields felt by bonds of interest, including subtraction of self-field contributions

## Prerequisites
* OpenBabel
* Gaussian 16
* GROMACS
* VMD
* Multiwfn
* Sobtop
* Python


### Required Python modules
* os
* subprocess
* Path
* defaultdict
* glob
* argparse
* math
* numpy
* pandas
* matplotlib
* re
* shutil
* scipy
* sys
* csv


## Organization
The main code ("```batch_auto_md_v#.sh```") should be placed in a folder with the solute molecule to be parametrized and intended solvents in separate folders. The structure should be:

```
Folder/
├── batch_auto_md_v#.sh
├── solute.cdx
└── Solvent1/
    └── SL1_100.gro
    └── SL1_100.top
└── Solvent2/
    └── SL2_100.gro
    └── SL2_100.top
└── 60Solvent3-40Solvent4/
    └── SL3_060.gro
    └── SL3_060.top
    └── SL4_040.gro
    └── SL4_040.top
```


## Usage
The main script ("```batch_auto_md_v#.sh```") needs to be able to reference the Starting_template folder to prepare the final folder structure, which will look like:

```
Folder/
├── batch_auto_md_v#.sh
└── 00_inputs/
    └── solute.cdx
└── 01_QM/
    └── (many files, .chk, .chg, .fchk, .mol, .mol2, etc.)
└── 02_resp/
    └── sobtop/
    └── all.resp
    └── VAC.resp
    └── water.resp
    └── Molecule_SP_VAC.chg
    └── Molecule_SP_water.chg
    └── 02-resp.py
└── 03_solute_params/
    └── all.resp
    └── Molecule-solv.itp
    └── Molecule.gro
    └── Molecule.itp
    └── Molecule.top
    └── sys.itp
    └── temp.top
    └── 03-extract_top_itp_sections.py
    └── 04-write_charges.py
└── Solvent1/
    └── 04_solvent_md/
    └── 05_production_md/
    └── 06_final_results/
    └── SL1_100.gro
    └── SL1_100.top
└── Solvent2/
    └── 04_solvent_md/
    └── 05_production_md/
    └── 06_final_results/
    └── SL2_100.gro
    └── SL2_100.top
└── 60Solvent3-40Solvent4/
    └── 04_solvent_md/
    └── 05_production_md/
    └── 06_final_results/
    └── SL3_060.gro
    └── SL3_060.top
    └── SL4_040.gro
    └── SL4_040.top
```

Final results files are copied to ```06_final_results/``` at the end of each solvent run, and electric field calculations are also performed and saved there.

# How to cite
If you found any of these functions useful, please consider citing this code as:

1. PA Kocheril et al. (in preparation).
