#bin/bash
module load gcc
module load openmpi
mpic++ -O3 main.cpp -o main
sbatch serial.sh
sbatch submit1.sh
sbatch submit4.sh
sbatch submit9.sh
sbatch submit16.sh
sbatch submit25.sh
sbatch submit36.sh
