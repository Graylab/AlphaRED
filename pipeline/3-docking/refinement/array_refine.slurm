#!/bin/bash -l
#----------------------------------------------------
# SLURM job script to run RosettaDock MPI application 
# with an array job to refine all structures.
#
# usage: from within the directory where your flag files are,
# submit to cluster using sbatch 
#----------------------------------------------------
#SBATCH --job-name=HR-dock
#SBATCH --time=6:00:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error logs/docking.%j.err       # Name of stdout output file (%j expands to jobId), remember to make outerr directory
#SBATCH --array=1-8
#SBATCH --output=logs/array_job_%A_%a.out

ml gcc

# job description
ROSETTABIN=$WORK/Rosetta/main/source/bin
ROSETTAEXE=docking_protocol
COMPILER=mpi.linuxgccrelease
EXE=$ROSETTABIN/$ROSETTAEXE.$COMPILER
echo Starting MPI job running $EXE
args=("decoys_P_0001_traj.out" "decoys_P_0002_traj.out" "decoys_P_0003_traj.out" "decoys_P_0004_traj.out" "decoys_P_0005_traj.out" "decoys_P_0006_traj.out" "decoys_P_0007_traj.out" "decoys_P_0008_traj.out")

# Retrieve the argument for the current array index
current_arg=${args[$SLURM_ARRAY_TASK_ID - 1]}

# Print the current argument
echo "Running job with argument: $current_arg"
date
time mpirun $EXE -in:file:silent $current_arg @local_refine_flags