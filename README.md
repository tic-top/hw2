# EECS87 Parallel computing HW2

## Usage
run all tests on greatlake.
```
bash run.sh 
```

Example, `mpirun -np 36 --bind-to core:overload-allowed ./main 2000 500 0 1 > parallel-small-36.txt` run a (2000,500) with 36 core on without interdemiate detail with parallel mode and output to  parallel-small-36.txt.

```
mpirun -np [number of core] --bind-to core:overload-allowed ./main [m] [n] [verbose(1) or not(0)] [parallel(1) or serial(0)] > [output path]
```

## Methos
I crop the matrix into $p$ blocks and each part has size $(\frac{m}{sqrt(p)}, \frac{n}{sqrt(p)})$. Except for the blocks on the boundary, the inner blocks need to send its left col, right col, bottom row, up row and 4 corners to the neibours.
## Result
I have test the model with more different numbers of processor(more than the problem mentiond.)

### Small (2000,500)
|    |   Processors |   Time Used |         Sum |   Sum of Square |   Speedup |   Efficiency |
|---:|-------------:|------------:|------------:|----------------:|----------:|-------------:|
|  0 |            1 |    6.1323   | 3.01222e+07 |     6.85949e+09 |  1        |     1        |
|  1 |            1 |    6.25376  | 3.01222e+07 |     6.85949e+09 |  0.980578 |     0.980578 |
|  2 |            4 |    1.57567  | 3.01222e+07 |     6.85949e+09 |  3.89187  |     0.972967 |
|  3 |           16 |    0.43552  | 3.01222e+07 |     6.85949e+09 | 14.0804   |     0.880026 |
|  4 |           36 |    0.191771 | 3.01221e+07 |     6.85949e+09 | 31.9772   |     0.888256 |

![image](small.png)

### Large (1000,4000)
|    |   Processors |   Time Used |         Sum |   Sum of Square |   Speedup |   Efficiency |
|---:|-------------:|------------:|------------:|----------------:|----------:|-------------:|
|  0 |            1 |   24.4224   | 1.19289e+08 |     1.41824e+10 |  1        |     1        |
|  1 |            1 |   24.7212   | 1.19289e+08 |     1.41824e+10 |  0.987913 |     0.987913 |
|  2 |            4 |    6.21868  | 1.19289e+08 |     1.41824e+10 |  3.92726  |     0.981816 |
|  3 |           16 |    1.60266  | 1.19289e+08 |     1.41824e+10 | 15.2387   |     0.952417 |
|  4 |           36 |    0.716405 | 1.19288e+08 |     1.41824e+10 | 34.0902   |     0.94695  |

![image](large.png)

### Scaling, Speedup and efficiency
From Strong Scaling perspective, we can see that the program is linear speedup but not perfect speedup since the efficiency is less than 1 as shown in the images above. I guess this maybe caused by the communication cost when the processor number goes up. So it's not perfect strong scaling.


From weak scaling perspective, we have scaling in terms of the size of the matrix, as well as in the number of processors. When we expand the matrix size 4 time, we compare the processor number parts like (1,4) (4,16) (9,36), we can see the running time is similar, so it's perfect weak scaling.