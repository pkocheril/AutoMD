#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=60:00:00   # walltime, maximum 168 hours
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "My First MD Simulation"   # job name
#SBATCH --mail-user=pkocheri@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


################## Setup ####################


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#MYHOME=`pwd`
#MYTMP="/central/scratch/$USER/scratch/$SLURM_JOBID"
#mkdir -p $MYTMP
#cp -r * $MYTMP/.
#cd $MYTMP
#echo "Start time at `date`"


################## Main script ####################

### 0. Rename starting file to Molecule.cdx

for file in *.cdx; do
    [ -f $file ] || continue
    mv $file Molecule.cdx
done

### 1. Copy starting template and move files

cp -r /central/groups/WeiLab/Phil/MD/Starting_template/* ./
mv *.* 00_inputs/
cd 00_inputs/

### 2. Build 3D geometry with OpenBabel

python 01-build3d_opt.py

#  - ensure LD_LIBRARY_PATH points to openbabel/lib64

#obabel Molecule.cdx -O ../01_convert/Molecule.smi
#obabel ../01_convert/Molecule.smi -O ../01_convert/Molecule_3d.mol
#obabel ../01_convert/Molecule_3d.mol -O ../02_geometry/Molecule_opt.mol

### 3. Run Gaussian calculations

cd ../03_gaussian/
python 02-gjf_maker.py
g16 Molecule_linked.gjf

### 4. Convert Gaussian output to .mol2 and make .fchk

for file in *.chk; do
    [ -f $file ] || continue
    formchk $file
    obabel ${file%.chk}.fchk -O ${file%.chk}.mol2
done

#python 03-gaussian_to_mol2.py

### 5. Run Multiwfn for RESP

for file in *.fchk; do
    [ -f $file ] || continue
    echo -e "7\n18\n1\ny\n" | Multiwfn $file
done

### 6. RESP

cp *SP*.chg ../04_resp/
cd ../04_resp/
python 04-resp.py

### 7. Parametrize with Sobtop

cp ../03_gaussian/Molecule_optfreq.fchk ../03_gaussian/Molecule_optfreq.mol2 ./sobtop/
cd sobtop/

# Run sobtop, using GAFF types and mSeminario for everything
echo -e "2\n\n1\n2\n2\nMolecule_optfreq.fchk\n\n\n0\n" | ./sobtop Molecule_optfreq.mol2

### 8. Move and modify molecule topology files

# -- noticed that OpenBabel's mol2 generation doesn't use the same numbering as GV

cp Molecule_optfreq.gro ../../05_solute_params/Molecule.gro
cp Molecule_optfreq.itp ../../05_solute_params/Molecule.itp
cp Molecule_optfreq.top ../../05_solute_params/Molecule.top
cd ../../05_solute_params/
python 05-extract_top_itp_sections.py
cp ../04_resp/all.resp ./
python 06-write_charges.py

### 9. Build a solvent box

cd ../06_solvent_md/
cp ../00_inputs/*.gro ../00_inputs/*.top ./
python 07-solvent_topology_edit.py

# Build solvent box with density specified in Solvent_properties.csv
python 08-build_box.py

# After this, you can visualize your solvent box with: vmd solv.gro

### 10. Solvent equilibration

# Manual version
#gmx grompp -f md/em.mdp -c solv.gro -p solv.top -o solv_em.tpr
#gmx mdrun -v -deffnm solv_em
#gmx grompp -f md/v.mdp -c solv_em.gro -p solv.top -o solv_v.tpr
#gmx mdrun -v -deffnm solv_v maxwarn 1
#gmx grompp -f md/p.mdp -c solv_v.gro -t solv_v.cpt -p solv.top -o solv_p.tpr
#gmx mdrun -v -deffnm solv_p
#gmx energy -f solv_v.edr -o temperature.xvg
#gmx energy -f solv_p.edr -o density.xvg

# Scripted version
bash full_auto_solvent.sh

# The above script runs:
#for file in solv*.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/em.mdp -c $file -p solv.top -o ${file%.gro}_em.tpr -maxwarn 0
#done
#for file in *_em.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -nt $(nproc)
#done
#for file in *_em.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/v.mdp -c $file -p solv.top -o ${file%_em.gro}_v.tpr -maxwarn 1
#done
#for file in *_v.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -cpo solv_v.cpt -nt $(nproc)
#done
#for file in *_v.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/p.mdp -c $file -t ${file%.gro}_v.cpt -p solv.top -o ${file%_v.gro}_p.tpr -maxwarn 0
#done
#for file in *_p.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -nt $(nproc)
#done

# Density and temperature equilibration
#for file in *_v.edr; do
#    [ -f $file ] || continue
#    echo -e "14\n\n" | gmx energy -f $file -o temperature.xvg
#done
#for file in *_p.edr; do
#    [ -f $file ] || continue
#    echo -e "22\n\n" | gmx energy -f $file -o density.xvg
#done
# After this, these equilibrations can be visualized with MATLAB (visualize_temperature_density.m)

### 11. Copy files for production MD

cd ../07_production_md/
cp ../05_solute_params/temp.top ../05_solute_params/Molecule.gro ../06_solvent_md/*_p.gro ./
# Rename file - weird naming previously
for file in *_p.gro; do
    [ -f $file ] || continue
    mv $file solv_p.gro
done
# mv solv.gro_p.gro solv_p.gro
mkdir itp
cp ../05_solute_params/sys.itp ../05_solute_params/Molecule-solv.itp ../06_solvent_md/*.itp ../06/solvent_md/solv.top ./itp
python 09-finalize_system_topology.py

# Make solvated box
gmx editconf -f Molecule.gro -o box0.gro -c -d 2.0 -bt cubic
gmx solvate -cp box0.gro -cs solv_p.gro -o MIX.gro -p topol.top

# Make 0q files
python 10-prepare-0q-topology.py


### 12. Production MD

# Scripted version
bash full_auto_prodmd.sh

# the above script runs:
#for file in *MIX.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/em.mdp -c $file -p topol.top -o ${file%.gro}_em.tpr -maxwarn 0
#done
#for file in *em.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
#done
#for file in *_em.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/v.mdp -c $file -p topol.top -o ${file%_em.gro}_v.tpr -maxwarn 0
#done
#for file in *v.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
#done
#for file in *v.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/p.mdp -c $file -t ${file%.gro}.cpt -p topol.top -o ${file%_v.gro}_p.tpr -maxwarn 0
#done
#for file in *p.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
#done
#for file in *p.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/md.mdp -c $file -t ${file%.gro}.cpt -p topol.top -o ${file%_p.gro}_md.tpr -maxwarn 0
#done
#for file in *p.gro; do
#    [ -f $file ] || continue
#    gmx grompp -f md/md.mdp -c $file -t ${file%.gro}.cpt -p topol_0q.top -o ${file%_p.gro}_0q.tpr -maxwarn 1
#done
#for file in *md.tpr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
#done
#for file in *md.trr; do
#    [ -f $file ] || continue
#    gmx mdrun -v -deffnm ${file%_md.trr}_0q -rerun ${file%.trr}.trr -nt ${nproc}
#done

# Trajectory conversion for periodic boundary conditions
gmx trjconv -s MIX_md.tpr -f MIX_md.xtc -o MIX_noPBC.xtc -pbc mol -center

# After this, you can use VMD to visualize the simulation trajectory

### 13. Post-run analysis

bash traj_extract.sh

# Move to final folder
cp *.xvg Molecule-solv.itp ../08_final_results/
cd ../08_final_results/

# Calculate electric fields
python 11-analyze_fields.py
# this will save a .png and .csv with the final field results

################## Cleanup ####################

#cp -r * $MYHOME
#cd $MYHOME
##rm -r $MYTMP
#echo "All ends at `date`" 

