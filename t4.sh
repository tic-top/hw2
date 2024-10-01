mpic++ main.cpp -o main
mpirun -np 4 --bind-to core:overload-allowed ./main 500 1000 0 > out4.txt
