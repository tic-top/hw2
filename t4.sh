mpic++ main.cpp -o main
mpirun -np 16 --bind-to core:overload-allowed ./main 2000 500 0 1 > out4.txt
