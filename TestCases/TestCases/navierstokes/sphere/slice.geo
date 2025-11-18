//+
SetFactory("OpenCASCADE");
// Set global mesh size
lc = 0.8;
//
Point(1) = {-1, 0, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {10, 0, 0, lc};
//+
Point(4) = {-10, 0, 0, lc};
//+
Point(5) = {0, 0, 0, lc};
//+
Circle(1) = {1, 5, 2};
//+
Circle(2) = {4, 5, 3};
//+
Line(3) = {4, 1};
//+
Line(4) = {2, 3};
//+
Curve Loop(1) = {2, 3, -1, 4};
//+
Plane Surface(1) = {1};
// Create a distance field and refine mesh around symmetry and sphere 
Field[1] = Distance;
Field[1].PointsList = {5};
Field[1].CurvesList = {1,3,4};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc / 5;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.15;
Field[2].DistMax = 0.5;

Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

//Mesh.MeshSizeExtendFromBoundary = 0;
//Mesh.MeshSizeFromPoints = 0;
//Mesh.MeshSizeFromCurvature = 0;

//Mesh.Algorithm = 5;
//+
//+
//Extrude {{1, 0, 0}, {1, 0, 0}, Pi/4} {
//  Surface{1}; Layers{5}; Recombine;
//}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/20} {
  Surface{1}; Layers{3}; Recombine;
}
//+
Physical Surface("symmetry_left", 11) = {4};
//+
Physical Surface("symmetry_right", 12) = {1};
//+
Physical Surface("farfield", 13) = {2};
//+
Physical Surface("wall", 14) = {3};
//+
Physical Volume("interior", 15) = {1};
