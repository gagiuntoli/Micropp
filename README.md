# Micropp

[![Build Status](https://travis-ci.org/GG1991/Micropp.svg?branch=master)](https://travis-ci.org/GG1991/Micropp)

Micropp is a 3-D Finite Element (FEM) code for solving solid
mechanics problems for modeling micro-scale effects in composite
materials.

The code solve the equilibrium equations at Representative Volume
Elements (RVE) to calculate average properties (stresses and
constitutive tensors) for multi-scale simulations:

<img src="./pics/mic_1.png" alt="drawing" width="300"/>

The code is mainly designed to be coupled with codes for
modeling with FEM the macro-scale, e.g. a wing of an aircraft.

<img src="./pics/coupling-micropp-macro.png" alt="drawing" width="300"/>

The code has been ported to GPUs to accelerate the calculation of the
micro-scale problems.

# Compilation 

```shell
cmake -B build \
    -DCMAKE_BUILD_TYPE=[Debug|Release] \
    -DENABLE_CUDA=[On|Off] \
    -DENABLE_OPENACC=[On|Off] \
    -DENABLE_OPENMP=[On|Off] \
    -DENABLE_TIMER=[On|Off]

cmake --build build
ctest --test-dir build
```

# Characteristics

1. Works with 3-D structured FE elements problems
2. OpenACC and CUDA (under development in `cuda` branch)
   acceleration support for GPUs
3. OpenMP support for multi-core CPUs
4. Solver: Diagonal Preconditioned Conjugate Gradients (DPCG)
5. More than 10 micro-structures patterns and 3 material laws
   (elastic, damage and plastic)
6. No external libraries are required
7. Native instrumentation to measure performance
8. C and Fortran Wrappers for coupling Micropp with external codes

# Performance CPU vs. GPUs

The peformance using the hybrid CPU/GPU execution in notably better
than CPU-only:

<img src="./pics/cpu-acc-cuda.png" alt="drawing" width="500"/>

For hybrid CPU/GPU the parallelization with CUDA has proven better
perfomance than OpenACC:

<img src="./pics/acc-cuda.png" alt="drawing" width="500"/>

Currently CUDA acceleration only works with some parts of the code and
has not been completely integrated.