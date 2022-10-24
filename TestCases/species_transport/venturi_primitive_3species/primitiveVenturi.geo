// ----------------------------------------------------------------------------------- //
// Tobias Kattmann, 09.09.2021, 2D Venturi Primitive for faster debugging
// ----------------------------------------------------------------------------------- //

// Evoque Meshing Algorithm?
Do_Meshing= 1; // 0=false, 1=true
// Write Mesh files in .su2 format
Write_mesh= 1; // 0=false, 1=true

// Geometric inputs
gas_inlet_diameter= 0.015;
air_inlet_diameter= 0.015;
gas_tube_length=0.045;
air_tube_length=0.045;
downstream_length= 0.09;

// Mesh sizing inputs
Nwall= 30; // Nodes for all spacings
gridsize= 0.1; // Later on not important as structured mesh is achieved

// ----------------------------------------------------------------------------------- //
// POINTS

// Starting in the origin, which is the most low-left point, and going clockwise.
// Gas inlet
Point(1) = {0, 0, 0, gridsize};
Point(2) = {0, gas_inlet_diameter, 0, gridsize};
//
Point(3) = {gas_tube_length, gas_inlet_diameter, 0, gridsize};
// Air inlet
Point(4) = {gas_tube_length, gas_inlet_diameter+air_tube_length, 0, gridsize};
Point(5) = {gas_tube_length+air_inlet_diameter, gas_inlet_diameter+air_tube_length, 0, gridsize};
//
Point(6) = {gas_tube_length+air_inlet_diameter, gas_inlet_diameter, 0, gridsize};
// outlet
Point(7) = {gas_tube_length+air_inlet_diameter+downstream_length, gas_inlet_diameter, 0, gridsize};
Point(8) = {gas_tube_length+air_inlet_diameter+downstream_length, 0, 0, gridsize};
//
Point(9) = {gas_tube_length+air_inlet_diameter, 0, 0, gridsize};
Point(10) = {gas_tube_length, 0, 0, gridsize};

// ----------------------------------------------------------------------------------- //
// LINES

// Gas inlet box, clockwise
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,10};
Line(4) = {10,1};

// air inlet box, clockwise
Line(5) = {3,4};
Line(6) = {4,5};
Line(7) = {5,6};
Line(8) = {6,3};

// downstream box, clockwise
Line(9) = {6,7};
Line(10) = {7,8};
Line(11) = {8,9};
Line(12) = {9,6};

// remaining lower middle box line
Line(13) = {9,10};

// ----------------------------------------------------------------------------------- //
// SURFACES (and Lineloops)
Curve Loop(1) = {1, 2, 3, 4}; Plane Surface(1) = {1}; // Gas inlet box
Curve Loop(2) = {5, 6, 7, 8}; Plane Surface(2) = {2}; // air inlet box
Curve Loop(3) = {9, 10, 11, 12}; Plane Surface(3) = {3}; // downstream box
Curve Loop(4) = {8, 3, -13, 12}; Plane Surface(4) = {4}; // middle box

// make structured mesh with transfinite Lines
// NOTE: The usage of Nwall and the progression has to be tuned again for any changes.
Transfinite Line{1,-3,12,-10} = Nwall Using Progression 0.9; // Spacing to the wall of the long tube, progression towards top wall
Transfinite Line{2,-4} = Nwall Using Progression 0.9; // Downstream spacing of gas inlet, progression towards air inlet

Transfinite Line{6,-8,-13} = Nwall Using Bump 0.1; // Spacing to the wall of the air inlet tube, progression towards side walls
Transfinite Line{-5,7} = Nwall Using Progression 0.9; // Downstream spacing of air inlet, progression towards gas inlet

Transfinite Line{-9,11} = Nwall Using Progression 0.9; // Downstream spacing of air inlet, progressio

// ----------------------------------------------------------------------------------- //
// PHYSICAL GROUPS

Physical Line("gas_inlet") = {1};
Physical Line("air_axial_inlet") = {6};
Physical Line("outlet") = {10};
Physical Line("axis") = {4,13,11};
Physical Line("wall") = {2,5,7,9};

Physical Surface("fluid") = {1,2,3,4};

// ----------------------------------------------------------------------------------- //
// Meshing
Transfinite Surface "*";
Recombine Surface "*";

If (Do_Meshing == 1)
    Mesh 1; Mesh 2;
EndIf

// ----------------------------------------------------------------------------------- //
// Write .su2 meshfile
If (Write_mesh == 1)

    Mesh.Format = 42; // .su2 mesh format, 
    Save "primitiveVenturi.su2";

EndIf

