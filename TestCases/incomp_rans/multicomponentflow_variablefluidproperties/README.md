# Test cases for multicomponent flow with composition dependent fluid properties  

## Introduction
Mixing laws for various fluid properties (viscosity, thermal conductivity, etc.) have been implemented in the SU2 software suite. This work is an extension of earlier work where equations for passive scalars and tabulated chemistry were implemented in the feature_flamelet branch.
The work in feature_multicomp deals with a transported scalar, i.e. there is a coupling between the scalar and flow solution. 
Both test cases concern the non-reacting mixing flow of methane and air through a planar (2D) venturi. The mass fraction of methane, which is the input to the mixing laws, is obtained by solving an additional transported scalar equation. 
 
The test cases concern:
- velocity inlets + pressure outlet
- pressure inlets + target mass flow rate 

Except for the boundary markers all other configuration options are the same.
The mesh file is also the same for both test cases. 

## Files for this test case
Below is a list of files for this test case and explanation for each file.
- SU2 repository
  - C6_pneumatic_venturi_planar_velocityinlets.cfg
    - Configuration file.
  - C6_pneumatic_venturi_planar_outlettargetmassflowrate.cfg
    - Configuration file.
- TestCases repository
  - C6_pneumatic_venturi_planar_mesh.cgns
    - Mesh file. The mesh was made with Fluent meshing. 

## How to use composition dependent fluid properties framework

### Config file
Enable composition dependent fluid properties using: 
- FLUID_MODEL = MIXTURE_FLUID_MODEL

For each fluid property specify an array of species’ fluid property e.g.:
- MOLECULAR_WEIGHT = 16.043, 28.965
  - Similar for specific heat, viscosity and thermal conductivity.

Define mixture diffusivity:
- DIFFUSIVITY_MODEL = UNITY_LEWIS
  - Mass diffusivity is computed based on unity Lewis assumption using mixture density, specific heat capacity (Cp) and thermal conductivity. 