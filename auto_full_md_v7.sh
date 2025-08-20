#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=24:00:00   # walltime, maximum 168 hours; at least 2 hours, probably 12-24
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks); at least 4, probably 16-20
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "Automated MD"   # job name
#SBATCH --mail-user=pkocheri@caltech.edu   # email address; update as needed
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


################## Setup ####################


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
MYHOME=`pwd`
MYTMP="/resnick/scratch/$USER/scratch/$SLURM_JOBID"
mkdir -p $MYTMP
cp -r * $MYTMP/.
cd $MYTMP
echo "Start time at `date`"


################## Main script ####################

##### Part 1: Parametrize molecule #####

### 1. Copy starting template and move files

cp -r /resnick/groups/WeiLab/Phil/MD/Automation/Starting_template/* ./
mv *.gro *.top *.cd* 00_inputs/
cd 01_QM/


### 2. Run Gaussian calculations

python 01-qm_prep.py --cores $(nproc) --mem $(nmem)


### 3. Run Multiwfn for RESP

for file in *.fchk; do
    [ -f $file ] || continue
    echo -e "7\n18\n1\ny\n" | Multiwfn $file
done


### 4. RESP

cp *SP*.chg ../02_resp/
cd ../02_resp/
python 02-resp.py

### 5. Parametrize with Sobtop

cp ../01_QM/Molecule_optfreq.fchk ../01_QM/Molecule_optfreq.mol2 ./sobtop/
cd sobtop/

# Run sobtop, using GAFF types and mSeminario for everything
echo -e "2\n\n1\n2\n2\nMolecule_optfreq.fchk\n\n\n0\n" | ./sobtop Molecule_optfreq.mol2


### 6. Move and modify molecule topology files

cp Molecule_optfreq.gro ../../03_solute_params/Molecule.gro
cp Molecule_optfreq.itp ../../03_solute_params/Molecule.itp
cp Molecule_optfreq.top ../../03_solute_params/Molecule.top
cd ../../03_solute_params/
python 03-extract_top_itp_sections.py
cp ../02_resp/all.resp ./
python 04-write_charges.py


##### Part 2: Loop over solvents #####

# Folders to ignore
F0="00_inputs/"
F1="01_QM/"
F2="02_resp/"
F3="03_solute_params/"
F4="04_solvent_md/"
F5="05_production_md/"
F6="06_final_results/"

for dir in */; do
	if [[ $dir == "$F0" ]] || \ 
	[[ $dir == "$F1" ]] ||
	[[ $dir == "$F2" ]] ||
	[[ $dir == "$F3" ]] ||
	[[ $dir == "$F4" ]] ||
	[[ $dir == "$F5" ]] ||
	[[ $dir == "$F6" ]] ; then
		continue # skip this folder
    fi
	echo "Processing directory: $dir"
	cd $dir
	
	
	### 7. Build a solvent box
	
	cp -r ../04_solvent_md/ ../05_production_md/ ../06_final_results/ ./
	cd 04_solvent_md/
	
	# Copy solvent .gro and .top files
	cp ../../00_inputs/*.gro ../../00_inputs/*.top ./
	python 05-solvent_topology_edit.py
	
	# Build solvent box with density specified in Solvent_properties.csv
	python 08-build_box.py
	# after this, you can visualize your solvent box with: vmd solv.gro
	
	
	### 8. Solvent equilibration
	
	bash full_auto.sh
	
	# Visualize density and temperature equilibrations
	python 07-plot_equilibration.py
	# this saves a figure called "solvent_equilibration.png" into 06_final_results
	
	
	### 9. Production MD
	
	cd ../05_production_md/
	cp ../03_solute_params/temp.top ../03_solute_params/Molecule.gro ../04_solvent_md/*_p.gro ./
	mkdir itp
	cp ../03_solute_params/sys.itp ../03_solute_params/Molecule-solv.itp ../04_solvent_md/*.itp ../04_solvent_md/solv.top ./itp/
	python 08-finalize_system_topology.py
	
	# Make solvated box
	gmx editconf -f Molecule.gro -o box0.gro -c -d 2.0 -bt cubic
	gmx solvate -cp box0.gro -cs solv_p.gro -o MIX.gro -p topol.top
	
	# Make 0q files
	python 09-prepare-0q-topology.py
	
	# Run MD
	bash full_auto_prodmd.sh
	
	# Trajectory conversion for periodic boundary conditions
	gmx trjconv -s MIX_md.tpr -f MIX_md.xtc -o MIX_noPBC.xtc -pbc mol -center
	# after this, you can use VMD to visualize the simulation trajectory
	
	
	### 10. Post-run analysis
	
	bash traj_extract.sh

	# Move to final folder
	cp *.xvg itp/Molecule-solv.itp ../06_final_results/
	cd ../06_final_results/

	# Calculate electric fields
	python 10-analyze_fields.py
	# this will save a .png and .csv with the final field histogram
	
	# Return to working directory for next loop iteration
	cd ../../
done


################## Cleanup ####################

## Remove extra template folders
#rm -r 04_solvent_md/ 05_production_md/ 06_final_results/

cp -r * $MYHOME
cd $MYHOME
#rm -r $MYTMP
echo "All ends at `date`" 

