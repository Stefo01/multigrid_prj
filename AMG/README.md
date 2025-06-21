# AMG: Algebraic Multigrid Solver

This repository contains a C++ implementation of an **Algebraic Multigrid (AMG) solver** for large sparse linear systems. AMG is a powerful iterative method for efficiently solving systems arising from discretized partial differential equations, especially those with large, sparse matrices.

## Features

- Classical AMG setup and solve phases
- Customizable coarsening and interpolation strategies
- FE implementation
- Sparse matrix storage using CSR format
- Gauss-Seidel smoothing
- Residual computation and monitoring
- Customizable mesh import
- Customizable boundary functions application
- VTU solution output

## Build Instructions

### Prerequisites

- C++17 or newer compiler (e.g., g++, clang++)
- CMake (recommended) or Make

### Build with CMake

```sh
AMG $ mkdir build
AMG $ cd build
build $ cmake ..
build $ make
```

## Usage

Run the solver with your input matrix and right-hand side:

```sh
build $ ./AMG
```

You can modify `main.cpp` to set up your problem, load matrices, and call the AMG solver.

You can run in parallel mode using:

```sh
build $ OMP_NUM_THREADS=1 ./AMG
```
with your desired number of threads.

## Main Classes

- **AMG**: Main class orchestrating the multigrid cycle.
- **CSRMatrix**: Compressed Sparse Row matrix storage and operations.
- **RestrictionOperator**: Handles coarsening and prolongation matrix construction.
- **Gauss_Seidel_iteration**: Smoother for each level.
- **FiniteElement**: Finite element discretization for the mesh
- **TriangularMesh**: Mesh managing

## Customization

- Adjust coarsening and interpolation strategies in `AMG.cpp` and `RestrictionOperator`.
- Change smoothing parameters in `apply_smoother_operator`.


