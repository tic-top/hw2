mpic++ main.cpp -o main
mpirun -np 1 --bind-to core:overload-allowed ./main 2000 500 0 0 > out1.txt
