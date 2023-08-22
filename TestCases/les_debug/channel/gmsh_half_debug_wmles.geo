//+
SetFactory("Built-in");
//+
SetFactory("OpenCASCADE");
//+ Set origin point
Point(1) = {0, 0, 0, 1.0};
//+ Set corner at x = 0, y = 0.1, z = 0
Point(2) = {0, 0.1, 0, 1.0};
//+ Define line 1 as line going from point 1 to point 2 (origin to delta)
Line(1) = {1, 2};
//+ Create a transfinite curve along line 1. The number of points should equal the number of desired layers plus 1. 
Transfinite Curve {1} = 6 Using Progression 1.5;
//+ Extrude the curves in the x-direction to get plane z_0
Extrude {2.5133, 0, 0} {
  Curve{1}; Layers {6}; Recombine;
  //Line{1,2}; Layers {136}; Recombine;
}
//+ Extrude plane z_0 in the z direction to get the volume
Extrude {0, 0, 0.9425} {
  Surface{1}; Layers {6}; Recombine;
}
//+ Name the surfaces to correspond to boundary conditions
//Physical Surface("x_0", 21) = {9, 5};
//+
//Physical Surface("y_0", 23) = {3};
//+
//Physical Surface("z_0", 22) = {2, 1};
//+
//Physical Surface("x_1", 25) = {10, 6};
//+
//Physical Surface("y_1", 24) = {8};
//+
//Physical Surface("z_1", 26) = {7, 11};
//+
Physical Surface("x_0", 13) = {4};
//+
Physical Surface("x_1", 14) = {5};
//+
Physical Surface("z_0", 15) = {1};
//+
Physical Surface("z_1", 16) = {6};
//+
Physical Surface("y_0", 17) = {2};
//+
Physical Surface("y_1", 18) = {3};
