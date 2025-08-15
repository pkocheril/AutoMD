#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=12:00:00   # walltime, maximum 168 hours; at least 2 hours, probably 12
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks); at least 4, probably 16
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "Automated MD"   # job name
#SBATCH --mail-user=pkocheri@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


################## Setup ####################


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#MYHOME=`pwd`
#MYTMP="/resnick/scratch/$USER/scratch/$SLURM_JOBID"
#mkdir -p $MYTMP
#cp -r * $MYTMP/.
#cd $MYTMP
#echo "Start time at `date`"


################## Main script ####################

cd Chloroform/
cp /resnick/groups/WeiLab/Phil/MD/auto_full_md_v6.sh ./
bash auto_full_md_v6.sh

cd ../DMSO/
cp /resnick/groups/WeiLab/Phil/MD/auto_full_md_v6.sh ./
bash auto_full_md_v6.sh

cd ../Hexane/
cp /resnick/groups/WeiLab/Phil/MD/auto_full_md_v6.sh ./
bash auto_full_md_v6.sh

cd ../THF/
cp /resnick/groups/WeiLab/Phil/MD/auto_full_md_v6.sh ./
bash auto_full_md_v6.sh

cd ../Toluene/
cp /resnick/groups/WeiLab/Phil/MD/auto_full_md_v6.sh ./
bash auto_full_md_v6.sh

cd ../Water/
cp /resnick/groups/WeiLab/Phil/MD/auto_full_md_v6.sh ./
bash auto_full_md_v6.sh


################## Cleanup ####################


#cp -r * $MYHOME
#cd $MYHOME
##rm -r $MYTMP
#echo "All ends at `date`" 

