mpic++ main.cpp -o main
mpirun -np 1 --bind-to core:overload-allowed ./main 300 300 0 > out.txt
