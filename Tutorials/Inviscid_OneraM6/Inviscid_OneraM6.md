Inviscid ONERA M6
=====

![ONERA M6 Cp](oneram6_cp.png)

## Goals

Upon completing this tutorial, the user will be familiar with performing a simulation of external, inviscid flow around a 3D geometry. The specific geometry chosen for the tutorial is the classic ONERA M6 wing. Consequently, the following capabilities of SU2 will be showcased in this tutorial:
- Steady, 3D Euler equations 
- Multigrid
- JST numerical scheme in space
- Euler implicit time integration
- Euler Wall, Symmetry, and Far-field boundary conditions
- Code parallelism (recommended)

We will also discuss the details for setting up 3D flow conditions and some of the multigrid options within the configuration file.

## Resources

The resources for this tutorial can be found in the TestCases/euler/oneram6/ directory. You will need the configuration file (inv_ONERAM6.cfg) and the mesh file (mesh_ONERAM6_inv_ffd.su2).

## Tutorial

The following tutorial will walk you through the steps required when solving for the flow around the ONERA M6 using SU2. The tutorial will also address procedures for both serial and parallel computations. It is assumed that you have already obtained and compiled the SU2_CFD code for a serial computation or both the SU2_CFD and SU2_SOL codes for a parallel computation. If you have yet to complete these requirements, please see the [[Download]] and [[Installation]] pages.

### Background

This test case is for the ONERA M6 wing in inviscid flow. The ONERA M6 wing was designed in 1972 by the ONERA Aerodynamics Department as an experimental geometry for studying three-dimensional, high Reynolds number flows with some complex flow phenomena (transonic shocks, shock-boundary layer interaction, separated flow). It has become a classic validation case for CFD codes due to the relatively simple geometry, complicated flow physics, and availability of experimental data. This test case will be performed in inviscid flow at a transonic Mach number.

### Problem Setup
This problem will solve the for the flow past the wing with these conditions:
- Freestream Pressure = 101325.0 N/m2
- Freestream Temperature = 288.15 K
- Freestream Mach number = 0.8395
- Angle of attack (AOA) = 3.06 deg

These transonic flow conditions will cause the typical "lambda" shock along the upper surface of the lifting wing.

### Mesh Description

The computational domain is a large parallelepiped with the wing half-span on one boundary in the x-z plane. The mesh provided in the tutorial resources directory listed above is a relatively coarse mesh provided that will complete in less time, but provide less accurate results. Users interested in obtaining more accurate results should use the finer mesh (mesh_ONERAM6_inv_FFD.su2) and associated config file (inv_oneram6_adv.cfg) provided in Testcases/optimization_euler/steady_oneram6/. The results shown in this tutorial use the finer mesh. 

The finer mesh consists of 582,752 tetrahedral elements and 108,396 nodes. Three boundary conditions are employed: Euler wall on the wing surface, a far-field characteristic-based condition on the far-field markers, and a symmetry boundary condition for the marker where the wing half-span is attached. The symmetry condition acts to mirror the flow about the x-z plane, reducing the complexity of the mesh and the computational cost. Images of the entire domain and the triangular elements on the wing surface are shown below.

![ONERA M6 Mesh](oneram6_mesh_bcs.png)
Figure (1): Far-field view of the computational mesh with boundary conditions.

![ONERA M6 Surface Mesh](oneram6_wing_mesh.png)
Figure (2): Close-up view of the unstructured mesh on the top surface of the ONERA M6 wing.

### Configuration File Options

Several of the key configuration file options for this simulation are highlighted here. The following describes how to set up 3D flow conditions in SU2:
```
% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%
%
% Mach number (non-dimensional, based on the free-stream values)
MACH_NUMBER= 0.8395
%
% Angle of attack (degrees)
AOA= 3.06
%
% Side-slip angle (degrees)
SIDESLIP_ANGLE= 0.0
%
% Free-stream pressure (101325.0 N/m^2 by default, only for Euler equations)
FREESTREAM_PRESSURE= 101325.0
%
% Free-stream temperature (288.15 K by default)
FREESTREAM_TEMPERATURE= 288.15
```

For an inviscid problem such as this, the flow conditions are completely defined by an input Mach number, flow direction, freestream pressure, and freestream temperature. The input Mach number is transonic at 0.8395. The freestream temperature and pressure have been set to standard sea level values for air at 101325.0 N/m2 and 288.15 K, respectively. The flow field will be initialized to these freestream values everywhere in the domain.

Lastly, it is very important to note the definition of the freestream flow direction in 3D. The default freestream direction (AOA = 0.0 degrees and SIDESLIP_ANGLE = 0.0 degrees) is along the positive x-axis without any components in the y- or z-directions. Referring to Figure (1), we see that AOA = 3.06 degrees will result in a non-zero freestream velocity in the positive z-direction. While zero for this problem, setting the SIDESLIP_ANGLE to a non-zero value would result in a non-zero velocity component in the y-direction. In 2D, the flow is in the x-y plane. While the default freestream direction is still along the positive x-axis, a non-zero AOA value for 2D problems will result in a non-zero freestream velocity in the y-direction. The SIDESLIP_ANGLE variable is unused in 2D.

In order to define reference values (for non-dimen. purposes):
```
% ---------------------- REFERENCE VALUE DEFINITION ---------------------------%
%
% Reference origin for moment computation
REF_ORIGIN_MOMENT_X = 0.25
REF_ORIGIN_MOMENT_Y = 0.00
REF_ORIGIN_MOMENT_Z = 0.00
%
% Reference length for pitching, rolling, and yawing non-dimensional moment
REF_LENGTH_MOMENT= 1.0
%
% Reference area for force coefficients (0 implies automatic calculation)
REF_AREA= 0
%
% Flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,
%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)
REF_DIMENSIONALIZATION= FREESTREAM_VEL_EQ_ONE
```

SU2 accepts arbitrary reference values for computing the force coefficients. A reference area can be supplied by the user for the calculation of force coefficients (e.g. a trapezoidal wing area) with the REF_AREA variable. If REF_AREA is set equal to zero, as for the ONERA M6, a reference area will be automatically calculated by summing all surface normal components in the positive z-direction on the monitored markers. For this ONERA M6 case SU2 performs a non-dimensional simulation (REF_DIMENSIONALIZATION= FREESTREAM_VEL_EQ_ONE). If you wish to perform a dimensional simulation you can pick the DIMENSIONAL option. For non-dimesionalization case FREESTREAM_PRESS_EQ_ONE the free-stream values at the farfield will be (pressure=1.0, density=1.0, temperature=1.0). For FREESTREAM_VEL_EQ_MACH the free-stream values at the farfield will be (velocity=MACH, density=1.0, temperature=1.0) and for FREESTREAM_VEL_EQ_ONE the free-stream values at the farfield will be (velocity=1.0, density=1.0, temperature=1.0).

Finally, we discuss some key multigrid options:
```
% -------------------------- MULTIGRID PARAMETERS -----------------------------%
%
% Multi-Grid Levels (0 = no multi-grid)
MGLEVEL= 3
%
% Multi-grid cycle (V_CYCLE, W_CYCLE, FULLMG_CYCLE)
MGCYCLE= W_CYCLE
```

SU2 contains an agglomeration multigrid algorithm for convergence acceleration technique where the original mesh is automatically agglomerated into a series of coarser representations, and calculations are performed on all mesh levels with each solver iteration in order to provide a better residual update. The user can set the number of multigrid levels using the MGLEVEL option. If this is set to zero, multigrid will be turned off, and only the original (fine) mesh will be used. An integer number of levels can be chosen. The ONERA M6 test case uses 3 levels of coarser meshes along with the original mesh for a total of 4 mesh levels. The type of cycle (V or W) can also be specified, and in general, while more computationally intensive, a W-cycle provides better convergence rates.

### Running SU2

Instructions for running this test case are given here for both serial and parallel computations. The computational mesh is rather large, so if possible, performing this case in parallel is recommended.

#### In Serial

The wing simulation is relatively large, but should still fit on a single-core machine. To run this test case, follow these steps at a terminal command line:
 1. Move to the directory containing the config file (inv_ONERAM6.cfg) and the mesh file (mesh_ONERAM6_inv_ffd.su2). Make sure that the SU2 tools were compiled, installed, and that their install location was added to your path.
 2. Run the executable by entering "SU2_CFD inv_ONERAM6.cfg" at the command line.
 3. SU2 will print residual updates with each iteration of the flow solver, and the simulation will terminate after reaching the specified convergence criteria.
 4. Files containing the results will be written upon exiting SU2. The flow solution can be visualized in ParaView (.vtk) or Tecplot (.dat for ASCII).

#### In Parallel

If SU2 has been built with parallel support, the recommended method for running a parallel simulation is through the use of the parallel_computation.py python script. This automatically handles the execution of SU2_CFD, and the writing of the solution files using SU2_SOL. Follow these steps to run the ONERA M6 case in parallel:
 1. Move to the directory containing the config file (inv_ONERAM6.cfg) and the mesh file (mesh_ONERAM6_inv_ffd.su2). Make sure that the SU2 tools were compiled, installed, and that their install location was added to your path.
 2. Run the python script by entering "parallel_computation.py -f inv_ONERAM6.cfg -n NP" at the command line with NP being the number of processors to be used for the simulation.
 3. SU2 will print residual updates with each iteration of the flow solver, and the simulation will terminate after reaching the specified convergence criteria.
 4. The python script will automatically call the SU2_SOL executable for merging the decomposed solution files from each processor into a single file. The files containing the results will be written upon exiting SU2. The flow solution can then be visualized in ParaView (.vtk) or Tecplot (.dat for ASCII).

### Results

Results are here given for the SU2 solution of inviscid flow over the ONERA M6 wing.

![ONERA M6 Cp](oneram6_cp.png)
Figure (3): Cp contours on the upper surface of the ONERA M6.

![ONERA M6 Mach](oneram6_mach.png)
Figure (4): Mach number contours on the upper surface of the ONERA M6 wing. Notice the "lambda" shock pattern typically seen on the upper surface.

![ONERA M6 Coefficients](oneram6_coefficients.png)
Figure (5): Convergence of the non-dimensional coefficients.

![ONERA M6 Convergence](oneram6_convergence.png)
Figure (6): Convergence of the density residual (speed up x20, iteration based).
