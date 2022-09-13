#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0051
#Outputs of the job
#SBATCH --output=OMPOutput.%j
#SBATCH --error=OMPError.%j
# Wall clock limit
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
## Specify the desired number of threads
#SBATCH --cpus-per-task=10
#SBATCH --exclusive
# run the process
./build-openmp-b/2d_Unsteady_OpenMP_B settings.coarse.in
