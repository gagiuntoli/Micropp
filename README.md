# _micropp_

Code to _localize_ strains and _homogenize_ stress in a Representative Volume Element (RVE) by using Finite Element Method (FEM).

# Characteristics

1. Works with structured grids 2D or 3D
2. Plastic non-linear material model for testing the memory storage and efficiency.
3. Supports boundary condition : uniform strains (Pure Dirichlet)
4. Runs sequentially.
5. Own ELL matrix routines with CG iterative solver (diagonal pre-conditioner).
6. Different kinds of micro-structures

# Main Characteristics

_micropp_ solves the FE problem on heterogeneous RVEs composed with more than one material and calculates the average properties of it. In the next figure a typical micro-structure is solved.

<img src="./pics/mic_1.png" alt="drawing" width="300"/>

_micropp_ is designed to be coupled with a macro-scale code in order to simulate multi-scale physical systems like an composite aircraft panel:

<img src="./pics/coupling-micropp-macro.png" alt="drawing" width="300"/>

`MicroPP` has been coupled with high-performance codes such as [Alya](http://bsccase02.bsc.es/alya) developed at the Barcelona Supercomputing center ([BSC](https://www.bsc.es/)) to performed **FE2** calculations. Also it was coupled with [MacroC](https://github.com/GG1991/macroc), a FE code that uses PETSc library on structured meshes. With this good performance was reach until 30720 processors on Marenostrum IV supercomputer.

<img src="./pics/scala.png" alt="drawing" width="350"/>

`MicroPP` has its own ELL matrix format routines optimized for the structured grid geometries that it has to manage. This allows to reach a really good performance in the assembly stage of the matrix. The relation between the assembly time and the solving time can be below than 1% depending on the problem size. The solving algorithm for the linear system of equations consists on a Conjugate Gradient algorithm with diagonal preconditioner.

Build steps with CMake:
-----------------------

1. Clone the repository 
2. cd cloned directory
3. mkdir build (can be also build+anything)
4. cd build
5. cmake .. (important the 2 points)
6. make

This will build the examples and the library in debug mode. CMake does not touch
the original sources and you can have many build directories with different
options if you want.

To build the optimized version:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

and the debug version:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Other possible options are: Debug, Release, RelWithDebInfo, MinSizeRel. Read CMake documentation for more information.

An option **TIMER** was added to insert time measure instrumentation for the execution. You can enable the option during cmake configuration time.

```bash
cmake -DTIMER=on ..
```

Option **CGDEBUG** was included to study the CG solver and see the convergence under different conditions.

```bash
cmake -DCGDEBUG=on ..
```

The new option is independent of Debug or release mode. But remember that any 


