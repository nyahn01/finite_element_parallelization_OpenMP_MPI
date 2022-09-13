****************************************************************************************************
OpenMP-PARALLEL VERSION OF 2D HEAT DIFFUSION EXPLICIT FEM SOLVER
****************************************************************************************************

The subject of this project will be to implement an OpenMP-parallel version of the solver of the heat equation
in 2D, which was implemented in Project 1 by means of the finite element method on an unstructured mesh. 

## Compilation

To compile the code you have to:

module load DEVELOP
module load intel/19.0          # The default is intel/19.0, but a later version is needed.
export CXX=`which icpc`
export C=`which icc`

mkdir build
cd build
cmake ..                        # Adjust accordingly the compiler flags in the CMakeLists.txt.
make

## Runnning

To run the code you have to:

cd test
sbatch < run.j                    # Adjust accordingly the run.j file, before submitting the job.

## Visualization

To visualize the result you have to:

module load GRAPHICS paraview
paraview

open disk-<coarse or fine or finest>.pvtu
