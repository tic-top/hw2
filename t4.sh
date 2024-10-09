mpic++ -O3 main.cpp -o main
mpirun -np 4 --bind-to core:overload-allowed ./main 2000 500 1 1 > out4.txt