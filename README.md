# MPI version for Tensor Train Decomposition
Parallelized C version of the [paper](https://arxiv.org/abs/1612.06505). 

# Structure of the code
- Data: contains input examples.
- One core: C version without MPI (no parallelization).
- Parallelization: C version with MPI (parallelized).

Each of the above directories contains:
- Algorithm1: Tensor Train (TT) inference/prediction.
- Algorithm2: Alternating Least Squares TT for training the model.

# Compiling the code
Remember to navigate to the algorithm path that you want to run, which must be one of the followings ```one_core/algorithm1```, ```one_core/algorithm2```, ```parallel/algorithm1```, ```parallel/algorithm2```:
```
cd /path/to/algorithm
make realclean 
make
```

# Input data

The executables may take one of the following input arguments (see examples in ```Data/```):
| Input | Description |
| ----------- | ----------- |
| poly_order.txt | File with number vandermonde vectors and dimensions |
| rank.txt | File with number of TT-ranks and each Grank  |
| X.txt | File with number of instances, dimension and all data variables |
| y.txt | File with number of instances and all data outputs variables (target/labels) |
| 48 | Matrix Blocksize for parallelizing <img src="https://latex.codecogs.com/gif.latex?C^t\cdot " /> <img src="https://latex.codecogs.com/gif.latex?C " /> (just for the  parallelized algorithm2)  |

Additionally, for algorithm1, the trained tensor weights (TT cores) and the input data are directly read from ```Data/algorithm1/G_*.txt``` and ```Data/X.txt```, respectively. Thus, you must save your TT cores in that path & format. 

# Running the code

### Algorithm1 - One core

```
./one_core/algorithm1/algorithm1 Data/algorithm1/poly_order.txt Data/algorithm1/rank.txt
```

### Algorithm2 - One core

```
./one_core/algorithm2/algorithm2 Data/algorithm2/poly_order.txt Data/algorithm2/rank.txt Data/algorithm2/X.txt Data/algorithm2/y.txt
```

### Algorithm1 - Parallel

```
mpiexec -n 4 ./parallel/algorithm1/algorithm1 Data/algorithm1/poly_order.txt Data/algorithm1/rank.txt Data/algorithm1/X.txt
```
where 4 specifies the number of cores.

### Algorithm2 - Parallel

```
mpiexec -n 4 ./parallel/algorithm2/algorithm2 Data/algorithm2/poly_order.txt Data/algorithm2/rank.txt Data/algorithm2/X.txt Data/algorithm2/y.txt 48
```
where 4 specifies the number of cores and 48 the matrix block size.
