// Gmsh .geo script
ri = 0.25;
ro = 0.5;
N = 40;

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {ri, 0.0, 0.0, 1.0};
Point(3) = {ro, 0.0, 0.0, 1.0};
Point(4) = {ro/Sqrt(2), ro/Sqrt(2), 0.0, 1.0};
Point(5) = {ri/Sqrt(2), ri/Sqrt(2), 0.0, 1.0};

Circle(1) = {2, 1, 5};
Circle(2) = {3, 1, 4};
Line(3) = {2, 3};
Line(4) = {5, 4};

Curve Loop(1) = {4, -2, -3, 1};
Surface(1) = {1};

Physical Curve("inlet", 5) = {1};
Physical Curve("outlet", 6) = {2};
Physical Curve("per1", 7) = {3};
Physical Curve("per2", 8) = {4};
Physical Surface("fluid", 9) = {1};

Transfinite Curve {1, 2} = N Using Progression 1;
Transfinite Curve {3, 4} = N Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

Mesh 2;
