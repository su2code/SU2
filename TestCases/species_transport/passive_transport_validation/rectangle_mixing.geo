// ----------------------------------------------------------------------------------- //
// Tobias Kattmann, 24.09.2021, 2D rectangle for mixing validation testcase
// ----------------------------------------------------------------------------------- //

// Evoque Meshing Algorithm?
Do_Meshing= 1; // 0=false, 1=true
// Write Mesh files in .su2 format
Write_mesh= 1; // 0=false, 1=true

// Geometric inputs
length= 10; // downstream direction
height= 5; // crossstream direction, full length
border= 1; // height of gas inlet, height-border= air inlet length

// Mesh sizing inputs. Note that the #cells+Prgression needs some fiddeling such that mesh size at the border fits.
Nl= 50; // Nodes in 'l'ength direction
Nhg= 25; // Nodes in 'h eight direction of the 'g'as side
Rhg= 0.9; // Progression of gas side [0,1], lower means more agressive
Nha= 40; // Nodes in 'h eight direction of the 'a'ir side
Rha= 0.9; // Progression of air side [0,1], lower means more agressive

gridsize= 0.1; // Later on not important as structured mesh is achieved

// ----------------------------------------------------------------------------------- //
// POINTS

// Starting in the origin, which is the most low-left point, and going clockwise.

Point(1) = {0, 0, 0, gridsize};
Point(2) = {0, border, 0, gridsize};
Point(3) = {0, height, 0, gridsize};
Point(4) = {length, height, 0, gridsize};
Point(5) = {length, border, 0, gridsize};
Point(6) = {length, 0, 0, gridsize};

// ----------------------------------------------------------------------------------- //
// LINES

// gas inlet 
Line(1) = {1,2};
// air inlet
Line(2) = {2,3};
// top sym
Line(3) = {3,4};
// air outlet
Line(4) = {4,5};
// gas outlet
Line(5) = {5,6};
// bottom sym
Line(6) = {6,1};
// species border
Line(7) = {2,5};

// ----------------------------------------------------------------------------------- //
// SURFACES (and Lineloops)
Curve Loop(1) = {6, 1, 7, 5};  Plane Surface(1) = {1}; // gas box
Curve Loop(2) = {3, 4, -7, 2}; Plane Surface(2) = {2}; // air box

// make structured mesh with transfinite Lines
// NOTE: The usage of Nwall and the progression has to be tuned again for any changes.
Transfinite Line{3,-6,7} = Nl ; // downstream direction, no progression
Transfinite Line{1,-5} = Nhg Using Progression Rhg; // gas side, progression towards border
Transfinite Line{-2,4} = Nha Using Progression Rha; // air side, progression towards border

// ----------------------------------------------------------------------------------- //
// PHYSICAL GROUPS

Physical Line("gas_inlet") = {1};
Physical Line("air_inlet") = {2};
Physical Line("outlet") = {4,5};
Physical Line("top") = {3};
Physical Line("bottom") = {6};

Physical Surface("fluid") = {1,2};

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
    Save "rectangle_mixing.su2";

EndIf
