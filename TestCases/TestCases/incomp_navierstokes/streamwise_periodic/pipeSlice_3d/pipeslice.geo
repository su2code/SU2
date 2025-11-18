//-------------------------------------------------------------------------------------//
// T. Kattmann, 13.05.2018, 3D Butterfly mesh in a circular pipe
// Create the mesh by calling this geo file with 'gmsh <this>.geo'.
//-------------------------------------------------------------------------------------//

// Evoque Meshing Algorithm?
Do_Meshing= 1; // 0=false, 1=true
// Write Mesh files in .su2 format
Write_mesh= 1; // 0=false, 1=true

//Geometric inputs, ch: channel, Pin center is origin
Radius= 0.5e-2; // Pipe Radius
InnerBox= Radius/2; // Distance to the inner Block of the butterfly mesh

//Mesh inputs
gridsize = 0.1; // unimportant once everything is structured

//ch_box 
Nbox = 30; // Inner Box points in x direction

Ncircu = 30; // Outer ring circu. points
Rcircu = 0.9; // Spacing towards wall

sqrtTwo = Cos(45*Pi/180);

//-------------------------------------------------------------------------------------//
//Points
// Inner Box
Point(1) = {-InnerBox, -InnerBox, 0, gridsize};
Point(2) = {-InnerBox, InnerBox, 0, gridsize};
Point(3) = {InnerBox, InnerBox, 0, gridsize};
Point(4) = {InnerBox, -InnerBox, 0, gridsize};

// Outer Ring
Point(5) = {-Radius*sqrtTwo, -Radius*sqrtTwo, 0, gridsize};
Point(6) = {-Radius*sqrtTwo, Radius*sqrtTwo, 0, gridsize};
Point(7) = {Radius*sqrtTwo, Radius*sqrtTwo, 0, gridsize};
Point(8) = {Radius*sqrtTwo, -Radius*sqrtTwo, 0, gridsize};

Point(9) = {0,0,0,gridsize}; // Helper Point for circles

//-------------------------------------------------------------------------------------//
//Lines
//Inner Box (clockwise)
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//Walls (clockwise)
Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};

//Connecting lines (outward facing)
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

//-------------------------------------------------------------------------------------//
//Lineloops and surfaces
// Inner Box (clockwise)
Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};

// Ring sections (clockwise starting at 9 o'clock)
Line Loop(2) = {5, -10, -1, 9}; Plane Surface(2) = {2};
Line Loop(3) = {10, 6, -11, -2}; Plane Surface(3) = {3};
Line Loop(4) = {-3, 11, 7, -12}; Plane Surface(4) = {4};
Line Loop(5) = {12, 8, -9, -4}; Plane Surface(5) = {5};

//make structured mesh with transfinite lines
//radial
Transfinite Line{1, 2, 3, 4, 5, 6, 7, 8} = Nbox;
//circumferential
Transfinite Line{9, 10, 11, 12} = Ncircu Using Progression Rcircu;

Transfinite Surface{1,2,3,4,5};
Recombine Surface{1,2,3,4,5};

//Extrude 1 mesh layer
Extrude {0, 0, 0.0005} {
    Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5};
    Layers{1};
    Recombine; 
}
Coherence;

//Physical groups made with GUI
Physical Surface("inlet") = {4, 1, 5, 3, 2};
Physical Surface("outlet") = {100, 122, 56, 78, 34};
Physical Surface("wall") = {69, 95, 113, 43};
Physical Volume("fluid") = {1, 2, 3, 4, 5};

// ----------------------------------------------------------------------------------- //
// Meshing
Transfinite Surface "*";
Recombine Surface "*";

If (Do_Meshing == 1)
    Mesh 1; Mesh 2; Mesh 3;
EndIf

// ----------------------------------------------------------------------------------- //
// Write .su2 meshfile
If (Write_mesh == 1)

    Mesh.Format = 42; // .su2 mesh format, 
    Save "pipe1cell3D.su2";

EndIf

