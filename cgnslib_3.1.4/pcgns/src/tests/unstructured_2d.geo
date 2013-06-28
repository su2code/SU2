Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};
Line Loop(5) = {1, 4, -2, -3};
Ruled Surface(6) = {5};
Recombine Surface {6};

Physical Surface('Internal') = {6};
Physical Line('South') = {1};
Physical Line('East') = {4};
Physical Line('North') = {2};
Physical Line('West') = {3};
