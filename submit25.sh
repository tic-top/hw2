#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)
# Set up batch job settings

#SBATCH --job-name=mpi_hw2
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --mem-per-cpu=1g
#SBATCH --time=00:05:00
#SBATCH --account=cse587f24_class
#SBATCH --partition=standard

mpirun -np 25 --bind-to core:overload-allowed ./main 2000 500 0 1 > parallel-small-16.txt
mpirun -np 25 --bind-to core:overload-allowed ./main 1000 4000 0 1 > parallel-large-16.txt