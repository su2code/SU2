//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.32, 0, 0, 1.0};
//+
Point(3) = {0.32, 0.0375, 0, 1.0};
//+
Point(4) = {0, 0.0375, 0, 1.0};
//+
Point(5) = {0, 0.0125, 0, 1.0};
//+
Point(6) = {-0.1, 0.0125, 0, 1.0};
//+
Point(7) = {-0.1, 0.0, 0, 1.0};
//+
Point(8) = {0.32, 0.0125, -0, 1.0};
//+
Line(1) = {7, 6};
//+
Line(2) = {1, 7};
//+
Line(3) = {1, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 8};
//+
Line(7) = {5, 8};
//+
Line(8) = {5, 4};
//+
Line(9) = {4, 3};
//+
Line(10) = {3, 8};
//+
Curve Loop(1) = {-2, 3, 4, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, -7, -3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {-7, 8, 9, 10};
//+
Plane Surface(3) = {3};
//+
Physical Curve("inlet", 11) = {1};
//+
Physical Curve("outlet", 12) = {6, 10};
//+
Physical Curve("wall_top", 13) = {9};
//+
Physical Curve("wall_side", 14) = {8};
//+
Physical Curve("wall_pipe", 15) = {4};
//+
Physical Curve("symmetry", 16) = {2, 5};
//+
Physical Surface("interior", 17) = {1, 2, 3};
//+
Transfinite Curve {2, 4} = 40 Using Progression 1.1;
//+
Transfinite Curve {5, 7, 9} = 100 Using Progression 1.04;
//+
Transfinite Curve {1, 3, 6} = 30 Using Progression 0.96;
//+
Transfinite Curve {8, 10} = 30 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Recombine Surface {1, 3, 2};
//+
Recombine Surface {1, 2, 3};
