## Micropp

Code to localize and average strain and stress over a micro structure.

# Characteristics

1. Structured grids 2D or 3D
2. Plasticity non-linear material model for testing the memory storage and efficiency.
3. Supports the three main kinds of boundary conditions: periodic, uniform strains and uniform stress.
4. Run on non-distributed architectures but can take advantage of multicore (OpenMP).
5. Can used external libraries (all that as we can) to solve the linear system of equations.

# Compile

Be sure of having the `g++` compiler and `make`, then run:
```bash
make
```
