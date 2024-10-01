# eecs587 hw2

```
module load gcc
module load openmpi
mpic++ main.cpp -o main
mpic++ hello.cpp -o hello_world
```

local

mpirun -np 1 ./main 500 500 0 > out.txt


mpirun -np 4 --bind-to core:overload-allowed ./hello_world > out.txt
