## Las Vegas Algorithm to solve ECDLP

This repository contains the implementation of the [Las Vegas Algorithm](http://issc.unipune.ac.in/attachments/lasVegas.pdf) in C++. Serial and Parallel implementations are available.
ECDLP_Parallel and ECDLP_Serial directories contains the parallel and serial implementations respectively.

### About the code.
Projective coordinates were used to implement Elliptic curve operations. Input for algorithm
is generated using SageMath. This algorithm is implemented for Binary as well as prime fields. 
The classes EC_GF2E and EC_ZZp implement the Las Vegas algorithm as their member function  which is named lasVegasECDLP.
The Las Vegas algorithm takes the following as input: Elliptic curve points P and Q, 
the order of P and a interger value called offset.
The serial implementation generates a single kernel and executes the 2-minor algorithm.


#### Parallel Implementation:
All the processors generate a kernel simultaneously. Each processor stores the kernel in a 
separate file. A directory named kernel stores these kernel files. The root processor reads the files 
from the kernel directory and computes partition data.
Partition data is the division of work among processors.
The root processor uses the Open MPI broadcast functionality to share the kernel along with the partition data among processors. 
Each processor upon receiving data executes the 2-minor algorithm.
The 2-minor algorithm searches for a 2x2 sub-matrix in the kernel such that the determinant is zero.
Searching for 2-minors is an embarrassingly parallel problem. 
Thus the speed up achieved by a parallel implementation over the serial implementation is substantial.


### Dependencies:
	1. C/C++ compiler
	2. OpemMPI
	3. Number Theory Library (NTL) and GNU-GMP
	4. Sage (genetate input)

### Compiling/Running Parallel Code.
	1. Parallel Code:- mpic++ *.cpp -lntl -lgmp -O3 -o ecdlp
	2. Running the Parallel Code : mpirun -n 2 ecdlp

### Compiling/Running Serial Code
	1. Serial Code:- g++ *.cpp -lntl -lgmp -O3 -o ecdlp
	2. Running Serial Code : ./ecdlp

### Compiling/Running Parallel code inside APM directory
	1. This directory is self contained. It consists all the code for algorithms to solve ECDLP.
	2. This directory also contains binary inputs for 20 thorough 80 bits. 
	3. Code here can be compiled using make command. An executable named lasVegas is generated.
	4. Run the code : mpirun -n 20 ./lasVegas