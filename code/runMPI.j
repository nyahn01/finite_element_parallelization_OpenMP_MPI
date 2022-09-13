#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0051
### Select the number of MPI PEs to use. 
#SBATCH --ntasks=8
#Outputs of the job
#SBATCH --output=MPIOutput.%j
#SBATCH --error=MPIError.%j
# Wall clock limit
#SBATCH --time=1:00:00
#SBATCH --exclusive
### run the process
### select a mesh
$MPIEXEC $FLAGS_MPI_BATCH ./build-mpi/2d_Unsteady_MPI settings.medium.in
