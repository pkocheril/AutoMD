#!/bin/bash
echo "Starting time at `date`"
echo "Running on hosts: $SLURM_NODELIST"
MYHOME=`pwd`
MYTMP="/central/scratch/$USER/scratch/$SLURM_JOBID"
mkdir -p $MYTMP
cp -r * $MYTMP/.
cd $MYTMP

run-your-job-script-here

cp * $MYHOME
cd $MYHOME
rm -r $MYTMP
echo "All ends at `date`" 
