# Miniflow

Miniflow is a mini-app testbed for 3D Navier-Stokes solution methods with various models of parallelism. The code and numerical methods are adapted from references [1] and [2].

## Build Instructions
Currently, miniflow has no automated build system. To build manually, use the following command:

### GCC
```g++ miniflow.cpp -O3 -march=native -std=c++20 -fopenmp -o miniflow```

### clang/LLVM
```clang++ miniflow.cpp -O3 -march=native -std=c++20 -fopenmp -o miniflow```

## Running Miniflow
Miniflow currently accepts TOML-formatted input files. Check out the ```examples/``` folder for example inputs.

## References
1. Griebel, M. <i>Numerical Simulation in Fluid Dynamics</i>. Society for Industrial and Applied Mathematics. 1998.
2. https://ins.uni-bonn.de/content/software-nast2d