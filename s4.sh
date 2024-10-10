#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)
# Set up batch job settings

#SBATCH --job-name=mpi_hw2
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=1g
#SBATCH --time=00:05:00
#SBATCH --account=cse587f24_class
#SBATCH --partition=standard

mpirun -np 4 --bind-to core:overload-allowed ./main 2000 500 1 1 > parallel-small-4.txt
mpirun -np 4 --bind-to core:overload-allowed ./main 1000 4000 1 1 > parallel-large-4.txt