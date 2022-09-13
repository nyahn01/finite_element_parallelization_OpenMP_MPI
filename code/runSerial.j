#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0051
#Outputs of the job
#SBATCH --output=SerialOutput.%j
#SBATCH --error=SerialError.%j
# Wall clock limit
#SBATCH --time=1:00:00
#SBATCH --exclusive

# run the process
./build-serial/2d_Unsteady_Serial settings.coarse.in
