// Structured wedge mesh for the Morse swirling jet (ERCOFTAC C26).
// The x-r plane is meshed at theta = -angle/2 and rotated about the jet axis (x) by the wedge angle, one cell thick.
// A slip wall at r = rmin replaces the axis, so all cells are hexahedra. Lengths are in nozzle diameters times D.
//
//   gmsh morse_wedge.geo -3 -format su2 -o morse_wedge.su2
//
// The mesh must be 3D (-3); a 2D mesh gives a file without volume elements.

D = 0.0246;        // nozzle diameter [m]
angle = 5;         // wedge angle [deg]
x0 = 0.5;          // inlet plane, x/D
xmax = 50;         // outlet plane, x/D
rmin = 0.02;       // slip wall replacing the axis, r/D
rjet = 0.5905;     // outer radius of the measured inlet profile, r/D
rfine = 4;         // end of the uniform radial spacing, r/D
rmax = 30;         // far-field radius, r/D
dr = 0.02;         // uniform radial spacing, r/D
dx0 = 0.02;        // axial spacing at the inlet, x/D
dxmax = 0.25;      // largest axial spacing, x/D
g = 1.03;          // growth ratio of stretched spacings

// Cells per segment: axial spacing grows from dx0 to dxmax, then stays uniform.
nx1 = Ceil(Log(dxmax / dx0) / Log(g));
xg = x0 + dx0 * (g^nx1 - 1) / (g - 1);
nx2 = Round((xmax - xg) / dxmax);
nr1 = Round((rjet - rmin) / dr);
nr2 = Round((rfine - rjet) / dr);
nr3 = Round(Log(1 + (rmax - rfine) * (g - 1) / dr) / Log(g));

xs[] = {x0, xg, xmax};
rs[] = {rmin, rjet, rfine, rmax};
nxs[] = {nx1, nx2};
gxs[] = {g, 1};
nrs[] = {nr1, nr2, nr3};
grs[] = {1, 1, g};
a = -0.5 * angle * Pi / 180;

For i In {0:2}
  For j In {0:3}
    p = newp;
    Point(p) = {xs[i] * D, rs[j] * D * Cos(a), rs[j] * D * Sin(a)};
    P[i * 4 + j] = p;
  EndFor
EndFor

// Lines of constant r (H) and of constant x (V), all oriented towards larger x or r.
For i In {0:1}
  For j In {0:3}
    l = newl;
    Line(l) = {P[i * 4 + j], P[(i + 1) * 4 + j]};
    Transfinite Curve{l} = nxs[i] + 1 Using Progression gxs[i];
    H[i * 4 + j] = l;
  EndFor
EndFor
For i In {0:2}
  For j In {0:2}
    l = newl;
    Line(l) = {P[i * 4 + j], P[i * 4 + j + 1]};
    Transfinite Curve{l} = nrs[j] + 1 Using Progression grs[j];
    V[i * 3 + j] = l;
  EndFor
EndFor

For i In {0:1}
  For j In {0:2}
    cl = newcl;
    Curve Loop(cl) = {H[i * 4 + j], V[(i + 1) * 3 + j], -H[i * 4 + j + 1], -V[i * 3 + j]};
    s = news;
    Plane Surface(s) = {cl};
    Transfinite Surface{s};
    Recombine Surface{s};
    S[i * 3 + j] = s;
  EndFor
EndFor

ext[] = Extrude {{1, 0, 0}, {0, 0, 0}, angle * Pi / 180} { Surface{S[]}; Layers{1}; Recombine; };

// Extrude returns, per source surface, the rotated copy, the volume and the four side surfaces.
For k In {0:5}
  top[k] = ext[6 * k];
  vol[k] = ext[6 * k + 1];
EndFor

e = 1e-4 * D;
big = 2 * rmax * D;
inlet[] = Surface In BoundingBox{x0 * D - e, -big, -big, x0 * D + e, big, big};
jet[] = Surface In BoundingBox{x0 * D - e, -big, -big, x0 * D + e, rjet * D + e, big};
outlet[] = Surface In BoundingBox{xmax * D - e, -big, -big, xmax * D + e, big, big};
axis[] = Surface In BoundingBox{x0 * D - e, -big, -big, xmax * D + e, rmin * D + e, big};
side[] = Surface In BoundingBox{x0 * D - e, rmax * D * Cos(a) - e, -big, xmax * D + e, rmax * D + e, big};
open[] = inlet[];
open[] -= jet[];

Physical Surface("inlet_jet") = {jet[]};
Physical Surface("farfield") = {open[], side[]};
Physical Surface("outlet") = {outlet[]};
Physical Surface("axis") = {axis[]};
Physical Surface("periodic_1") = {S[]};
Physical Surface("periodic_2") = {top[]};
Physical Volume("fluid") = {vol[]};
