#!/bin/bash
#SBATCH --job-name=abaqus-easypbc4inp
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --partition=operation
#SBATCH -t 20-12:8:8
#SBATCH --mem=800GB
#SBATCH -o %x_%j.log
### Configure environment variables, need to unset SLURM's Global Task ID for ABAQUS's PlatformMPI to work
### Create ABAQUS environment file for current job, you can set/add your own options (Python syntax)
module purge
module add abaqus/2022
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
ulimit -s unlimited

abq2022 cae noGUI=abaqus python full_pipeline_v1.2.py

