#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)
# Set up batch job settings

#SBATCH --job-name=mpi_hw2
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=1g
#SBATCH --time=00:05:00
#SBATCH --account=linmacse2
#SBATCH --partition=standard

# module load gcc
# module load openmpi
# mpic++ main.cpp -o main
mpirun -np 4 --bind-to core:overload-allowed ./main 300 1000 0 > test.txt
