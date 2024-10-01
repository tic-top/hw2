mpic++ main.cpp -o main
mpirun -np 4 --bind-to core:overload-allowed ./main 200 50 1 > out.txt
