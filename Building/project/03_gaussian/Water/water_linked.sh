#!/bin/bash

# Submit this script with: sbatch water_linked.sh

#SBATCH --job-name=water_linked
#SBATCH --time=167:59:59
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-user=pkocheri@caltech.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

g16 water_linked.gjf
