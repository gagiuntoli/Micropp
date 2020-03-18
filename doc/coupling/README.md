# README

This folder contains Alya inputs for simulating FE2 multi-scale
problems and use the results for the paper. For executing these
simulations Alya should be compiled and linked together with Micropp
code in order to solve FE2.

# Compilation

    tar -xvzf micropp.tar.gz
    tar -xvzf alyafe2.tar.gz

    cd micropp
    mkdir build-cpu
    mkdir build-gpu
    
    cd build-cpu
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    ctest (check everything is fine 100% of the cases should pass)

    cd ..

    cd build-gpu
    cmake .. -DCMAKE_BUILD_TYPE=Release -DOPENACC=On
    make
    ctest (check everything is fine 100% of the cases should pass)

    cd alyafe2/Thirdparties/micropp and create the symbolic links:

    cd build-cpu
    ln -s ../../../../micropp/build-cpu/libmicropp.a
    ln -s ../../../../micropp/build-cpu/libmicropp.mod
    ln -s ../../../../micropp/build-cpu/libmaterials.mod

    cd ..

    cd build-gpu
    ln -s ../../../../micropp/build-gpu/libmicropp.a
    ln -s ../../../../micropp/build-gpu/libmicropp.mod
    ln -s ../../../../micropp/build-gpu/libmaterials.mod

    cd ../../../.. 

    cd Executables/cpu
    ./configure solidz parall
    make -j4 MICROPP=1

    cd ../gpu
    ./configure solidz parall
    make -j4 MICROPP=1


# Basic FE2 Simulation

The geometry used in the macro-scale is a prismatic piece with a
fracture in order to concentrate tension on a specific region of the
volume and to have a localized area where the micro-scale enters in
the non-linear zone. On the other hand the micro-scale is a cube with
porouses distributed randomly.

All the elements of the macro-scale mesh are coupled with 8
micro-scale problems. The macro-scale transfers the strain tensor to
the micro-scale, this last uses the strain to compute the equilibrium
in the cubic cell and returns the stress tensor and the tangent
constitutive tensors to the macro-scale.

The problem becomes non-linear if at least one of the micro-scale
problems enter in the non-linear region. This happens if at least one
of the integration points of the micro-scale modeled with a non-linear
material passes a strain or stress limit.

In a tipical problem the macro-scale can contain certain regions with:

 * No-Coupling coupling, i.e., materials will remain linear in the
   entire simulation (computationally cheap).

 * One-Way Coupling: computes a FE problem for getting the stress.

 * Full Coupling: solves 7 FE problem for getting the stress and the
   tangent constitutive tensors (computationally expensive).

# Alya Input

Alya inputs specified how the macro-scale is executed: time step,
solver used and other important parameters used in solid mechanics.
case the micro-scale is a constitutive law.

In the inputs of this case the file that determines the caracteristics
of the micro-structure and its simulation parameters are the *.sld.dat
files (one per simulation).

The parameters of the micro-struture should be set in first material.
Note that the current implementation does not allow to set a wide
variety of micro-structure on the same simulation.

The following format is currently adopted:

    CONSTITUTIVE_MODEL   MICNC   <nx>  <ny>  <nz>   <micro-type> \
    <p1> <p2> <p3> <p4> <mat-1> <mat-2> <mat-3> \
    <E1>  <poisson1>  <Sy1>  <K1>  <H1> \
    <E2>  <poisson2>  <Sy2>  <K2>  <H2> \
    <E3>  <poisson3>  <Sy3>  <K3>  <H3> \
    <subiterations>

where (all parameters refer to the micro-scale):

 * `<nx,y,z>`: are the number of nodes in each of the directions of the
cell.

 * `<micro-type>`: a number between 0 and 10 defining the kind of
micro-structure pattern.

 * `<p1,2,3,4>`: real numbers to modify the micro-structure, i.e. the
carbon fiber radious.

 * `<mat-1,2,3>`: material kinds used for each constituents in the
micro-structure. Elastic (0), Plastic (1), Damage (2).

 * `<E>  <poisson>  <Sy>  <K>  <H>`: Are the constants that define
the material model.

Example:

    PROPERTIES
        $---------------------------------------
        MATERIAL           = 1
        DENSITY            = 1.6E-9
        CONSTITUTIVE_MODEL   MICNC \
        3  3  3   2   1.0  1.0  1.0  0.2  0  0  0 \
        3.0e7  0.25  0.0e0 4.0e3  4.0e3 \
        3.0e7  0.25  0.0e0 4.0e3  4.0e3 \
        3.0e7  0.25  0.0e0 4.0e3  4.0e3 \
        1
        $---------------------------------------
        MATERIAL           = 2
        DENSITY            = 1.6E-9
        CONSTITUTIVE_MODEL   MICOW
        $---------------------------------------
    END_PROPERTIES

This mean that there are two material in the mesh. The one with `ID =
1` will have the properties of `MATERIAL = 1` and the other with `ID =
2` will have the properties of `MATERIAL = 2`.

For `MATERIAL = 2` the model `MICNC` means Micropp No-Coupling this
means that Micropp assumes that the micro-structure will remain in the
linear zone during the entire simulation. In `MATERIAL = 2` the
Micropp model is `MICOW` meaning Micropp One-Way, this model only
computes the stress tensor applying 1 FE computation.


# Parallel Execution

Alya uses domain partition methods to solve the macro-scale in
parallel. Each MPI computes a non-overlapping region of the domain and
then the information in the boundaries between subdomain is transfered
across MPIs.

For FE2 the coupling of Alya with Micropp is made by instantiating one
Micropp object in each MPI of Alya. This instance of Micropp computes
all the micro-scale problem corresponding to the integration points in
that subdomain. Since the micro-scale computation can take lot of time
depending on the resolution of its mesh, it is important to have a
balanced distribution of micro-scale problems across MPIs in Alya.


Example with a macro-scale containing 2 hexahedron element (8 Gauss
points per elements).

   Execution: `mpirun -np 3 Alya input`

    MPI0 (Master)  -   MPI1   -    MPI2
                        |           |
                        |           |
                      Micropp      Micropp
   		        |           |
                       gp0         gp8
                       gp1         gp9
                       gp2         gp10
                       gp3         gp11
                       gp4         gp12
                       gp5         gp13
                       gp6         gp14
                       gp7         gp15


