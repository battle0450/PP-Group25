# PP-Group25
Parallel Programming Final Project - Group25. <br>

How to Compile:
`make`

Exeucting PSO_serial:
```C
./PSO_serial N D ITER
// N = Num of particles
// D = Dimensions
// ITER = Num of Iterations
```

Executing PSO_pthread:
```C
./PSO_pthread N D ITER THREADS
// N = Num of particles
// D = Dimensions
// ITER = Num of Iterations
// THREADS = Num of threads
```

Executing PSO_mpi:
```C
mpirun -np P --npernode 1 --hostfile hosts ./mpi N D ITER
// P = Num of processes
// N = Num of particles
// D = Dimensions
// ITER = Num of Iterations
```

Executing PSO_CUDA:
```C
./PSO
//Modify N D ITER in kernel.h
```
