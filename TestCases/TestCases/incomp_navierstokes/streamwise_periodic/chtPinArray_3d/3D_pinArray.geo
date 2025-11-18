// ------------------------------------------------------------------------- //
// T. Kattmann, 18.06.2019, 3D 2 Zone mesh
// Create the mesh by calling this geo file with 'gmsh <this>.geo'.
// For multizone mesh the zonal meshes have to be created using the first 
// option 'Which_Mesh_Part' below and have to be married appropriatley.
// ------------------------------------------------------------------------- //

// Which domain part should be handled
Which_Mesh_Part= 1; // 0=all, 1=Fluid, 2=Solid, 3=InterfaceOnly
// Evoque Meshing Algorithm?
Do_Meshing= 1; // 0=false, 1=true
// Write Mesh files in .su2 format
Write_mesh= 1; // 0=false, 1=true
// Mesh Resolution
Mesh_Resolution= 1; // 0=debugRes, 1=Res1, 2=Res2
// Have conical pins?
ConicalPins= 1; // 0=No, i.e. cylindrical, 1=Yes 

// Free parameters
scale_factor= 1e-3; // scales Point positions from [mm] to [m] with 1e-3
dist= 6.44 * scale_factor; // distance between pin midpoints, each pin has 6 surrounding pins, i.e. 60 deg between each
upper_h= 10.0 * scale_factor; // height (z-dir) of fluid domain/pins
lower_h= 5.0 * scale_factor; // height (z-dir) of solid domain/pins

If (ConicalPins==0) // cylindrical
    r_pin_lower= 2.0 * scale_factor; // lower pin radius
    r_pin_upper= 2.0 * scale_factor; // upper pin radius
ElseIf(ConicalPins==1) // conical
    r_pin_lower= 2.3 * scale_factor; // lower pin radius
    r_pin_upper= 0.56 * scale_factor; // upper pin radius
EndIf

// Dependent parameters
rad2deg= Pi/180; // conversion factor as gmsh Cos/Sin functions take radian values
length= 2 * Cos(30*rad2deg)*dist; // length (in x-dir)
width= Sin(30*rad2deg)*dist; // width (in y-dir)

Printf("===================================");
Printf("Free parameters:");
Printf("-> distance between pins: %g", dist);
Printf("-> lower pin radius: %g", r_pin_lower);
Printf("-> upper pin radius: %g", r_pin_upper);
Printf("Dependent parameters");
Printf("-> length: %g", length);
Printf("-> width: %g", width);
Printf("Fixed parameters");
Printf("-> pin height: %g", upper_h);
Printf("-> solid thickness: %g", lower_h);
Printf("===================================");

// Mesh inputs
gs = 0.5 *scale_factor; // gridsize

If(Mesh_Resolution==0) // debugRes
    // interface meshing parameteres. Also sufficient for fluid domain meshing.
    N_x_flow= 10; // #gridpoints in flow x-direction on a patch. Also N_x_flow/2 on smaller patches employed.

    N_y_flow = 10; // #gridpoints normal to pin surface, y-direction
    R_y_flow= 1.08; // Progression normal to pin surface

    N_z_flow= 20; // #gridpoints in height z-direction
    R_z_flow= 0.05; // Bump in height as top and bottom are walls

    // Additional meshing parameters for solid domain
    InnerRadiusFactor= 0.6; // How much of the inner pin is unstructured mesh (0.9=mostly unstructured, 0.1= mostly structured). Requires 0 < value < 1. 
    N_y_innerPin= 5; // #gridpoints of the structured first part of the inner Pin in y-direction / normal to the pin
    R_y_innerPin= 0.93; // Progression towards interface
    N_z_solid= 10; // #points from bottom interface to heater surface
    R_z_solid= 1.18; // progression for N_z_solid

ElseIf(Mesh_Resolution==1) // Res1
    // interface meshing parameteres. Also sufficient for fluid domain meshing.
    N_x_flow= 20; // #gridpoints in flow x-direction on a patch. Also N_x_flow/2 on smaller patches employed.

    N_y_flow = 40; // #gridpoints normal to pin surface, y-direction
    R_y_flow= 1.08; // Progression normal to pin surface

    N_z_flow= 100; // #gridpoints in height z-direction
    R_z_flow= 0.05; // Bump in height as top and bottom are walls

    // Additional meshing parameters for solid domain
    InnerRadiusFactor= 0.6; // How much of the inner pin is unstructured mesh (0.9=mostly unstructured, 0.1= mostly structured). Requires 0 < value < 1. 
    N_y_innerPin= 30; // #gridpoints of the structured first part of the inner Pin in y-direction / normal to the pin
    R_y_innerPin= 0.93; // Progression towards interface
    N_z_solid= 20; // #points from bottom interface to heater surface
    R_z_solid= 1.18; // progression for N_z_solid

ElseIf(Mesh_Resolution==2) // Res2
    // interface meshing parameteres. Also sufficient for fluid domain meshing.
    N_x_flow= 30; // #gridpoints in flow x-direction on a patch. Also N_x_flow/2 on smaller patches employed.

    N_y_flow = 50; // #gridpoints normal to pin surface, y-direction
    R_y_flow= 1.08; // Progression normal to pin surface

    N_z_flow= 200; // #gridpoints in height z-direction
    R_z_flow= 0.05; // Bump in height as top and bottom are walls

    // Additional meshing parameters for solid domain
    InnerRadiusFactor= 0.6; // How much of the inner pin is unstructured mesh (0.9=mostly unstructured, 0.1= mostly structured). Requires 0 < value < 1. 
    N_y_innerPin= 40; // #gridpoints of the structured first part of the inner Pin in y-direction / normal to the pin
    R_y_innerPin= 0.93; // Progression towards interface
    N_z_solid= 30; // #points from bottom interface to heater surface
    R_z_solid= 1.18; // progression for N_z_solid

EndIf

// Feasability checks
If (r_pin_lower >= width || 
    r_pin_upper >= width ||
    r_pin_lower <= 0 ||
    r_pin_upper <= 0)
    Printf("Aborting! Bad Inputs");
    Abort;
EndIf

// ------------------------------------------------------------------------- //
// Point distribution rules:
// Line/curve orientation rules (in that order): 
//   1. always pointing in positive z-direction
//   2. in x-y-plane, always pointing in positive x-direction
//   3. in y-direction, always pointing in positive y-direction
// Surface normal orientation rules: none

// ------------------------------------------------------------------------- //
// CHT Interface, complete description as it is part of fluid and solid
// Id's starting with in the range (1-99)

// Points
// Lower Pin1
Point(10) = {0, width, 0, gs}; // lower pin1 midpoint
Point(11) = {0, width-r_pin_lower, 0, gs}; // lower pin1 on inlet
Point(12) = {Sin(30*rad2deg)*r_pin_lower, width-Cos(30*rad2deg)*r_pin_lower, 0, gs}; // lower pin1 in between
Point(13) = {r_pin_lower, width, 0, gs}; // lower pin1 on sym
Circle(10) = {11,10,12}; // lower pin1 smaller first part
Circle(11) = {12,10,13}; // lower pin1 larger second part

// Lower Pin2
Point(20) = {0.5*length, 0, 0, gs}; // pin midpoint
Point(21) = {0.5*length - r_pin_lower, 0, 0, gs}; // lower small x
Point(22) = {length/2 - Sin(30*rad2deg)*r_pin_lower, Cos(30*rad2deg)*r_pin_lower, 0, gs}; // small intermediate
Point(23) = {length/2 + Sin(30*rad2deg)*r_pin_lower, Cos(30*rad2deg)*r_pin_lower, 0, gs}; // large intermediate
Point(24) = {0.5*length + r_pin_lower, 0, 0, gs}; // lower large x
Circle(20) = {21,20,22}; // first segment
Circle(21) = {22,20,23}; // second segment
Circle(22) = {23,20,24}; // third segment

// lower Pin3
Point(30) = {length, width, 0, gs}; // midpoint
Point(31) = {length, width-r_pin_lower, 0, gs}; // on outlet
Point(32) = {length-Sin(30*rad2deg)*r_pin_lower, width-Cos(30*rad2deg)*r_pin_lower,0, gs};
Point(33) = {length - r_pin_lower, width, 0, gs}; // on sym
Circle(30) = {31,30,32}; // first segment
Circle(31) = {32,30,33}; // second segment

// lower additional structured mesh points
Point(40) = {length/4 + Tan(30*rad2deg)*width/2, width, 0, gs}; // first half, large y
Point(41) = {length/4 - Tan(30*rad2deg)*width/2, 0, 0, gs}; // first half, small y
Point(42) = {length*3/4 - Tan(30*rad2deg)*width/2, width, 0, gs}; // second half, large y
Point(43) = {length*3/4 + Tan(30*rad2deg)*width/2, 0, 0, gs}; // second half, small y
Point(44) = {0, 0, 0, gs}; // corner point inlet
Point(45) = {length, 0, 0, gs}; // corner point outlet

// lower additional structured mesh lines
// outer boundary
Line(40) = {11, 44};
Line(41) = {44, 41};
Line(42) = {41, 21};
Line(43) = {43, 24};
Line(44) = {43, 45};
Line(45) = {45, 31};
Line(46) = {33, 42};
Line(47) = {42, 40};
Line(48) = {40, 13};
// inner lines
Line(49) = {41, 12};
Line(50) = {41, 40};
Line(51) = {22, 40};
Line(52) = {23, 42};
Line(53) = {42, 43};
Line(54) = {43, 32};

// line loops and surfaces on lower domain interface
Line Loop(10) = {40, 41, 49, -10}; Plane Surface(10) = {10};
Line Loop(11) = {-49, 50, 48, -11}; Plane Surface(11) = {11};
Line Loop(12) = {42, 20, 51, -50}; Plane Surface(12) = {12};
Line Loop(13) = {-51, 21, 52, 47}; Plane Surface(13) = {13};
Line Loop(14) = {53, 43, -22, 52}; Plane Surface(14) = {14};
Line Loop(15) = {53, 54, 31, 46}; Plane Surface(15) = {15};
Line Loop(16) = {44, 45, 30, -54}; Plane Surface(16) = {16};

// No progression in flow direction on the pin surface
Transfinite Line {11,50,20,47,21,22,53,31} = N_x_flow;
Transfinite Line {10,41,30,44} = N_x_flow/2;
// Progression normal to the pin surface
Transfinite Line {40, -49, -48, -42, 51, 52, -43, 46, -54, -45} = N_y_flow Using Progression R_y_flow; 

// Upper Pin1
Point(10+100) = {0, width, upper_h, gs}; // lower pin1 midpoint
Point(11+100) = {0, width-r_pin_upper, upper_h, gs}; // lower pin1 on inlet
Point(12+100) = {Sin(30*rad2deg)*r_pin_upper, width-Cos(30*rad2deg)*r_pin_upper, upper_h, gs}; // lower pin1 in between
Point(13+100) = {r_pin_upper, width, upper_h, gs}; // lower pin1 on sym
Circle(10+100) = {11+100,10+100,12+100}; // lower pin1 smaller first part
Circle(11+100) = {12+100,10+100,13+100}; // lower pin1 larger second part

// Upper Pin2
Point(20+100) = {0.5*length, 0, upper_h, gs}; // pin midpoint
Point(21+100) = {0.5*length - r_pin_upper, 0, upper_h, gs}; // lower small x
Point(22+100) = {length/2 - Sin(30*rad2deg)*r_pin_upper, Cos(30*rad2deg)*r_pin_upper, upper_h, gs}; // small intermediate
Point(23+100) = {length/2 + Sin(30*rad2deg)*r_pin_upper, Cos(30*rad2deg)*r_pin_upper, upper_h, gs}; // large intermediate
Point(24+100) = {0.5*length + r_pin_upper, 0, upper_h, gs}; // lower large x
Circle(20+100) = {21+100,20+100,22+100}; // first segment
Circle(21+100) = {22+100,20+100,23+100}; // second segment
Circle(22+100) = {23+100,20+100,24+100}; // third segment

// Upper Pin3
Point(30+100) = {length, width, upper_h, gs}; // midpoint
Point(31+100) = {length, width-r_pin_upper, upper_h, gs}; // on outlet
Point(32+100) = {length-Sin(30*rad2deg)*r_pin_upper, width-Cos(30*rad2deg)*r_pin_upper, upper_h, gs};
Point(33+100) = {length - r_pin_upper, width, upper_h, gs}; // on sym
Circle(30+100) = {31+100,30+100,32+100}; // first segment
Circle(31+100) = {32+100,30+100,33+100}; // second segment

// connection in height/z-direction
// Pin1
Line(61) = {11,11+100};
Line(62) = {12,12+100};
Line(63) = {13,13+100};

// Pin2
Line(71) = {21,21+100};
Line(72) = {22,22+100};
Line(73) = {23,23+100};
Line(74) = {24,24+100};

// Pin3
Line(81) = {31,31+100};
Line(82) = {32,32+100};
Line(83) = {33,33+100};

// Pin surfaces
Line Loop(20) = {61, 10+100, -62, -10}; Surface(20) = {20};
Line Loop(21) = {62, 11+100, -63, -11}; Surface(21) = {21};
Line Loop(22) = {71, 20+100, -72, -20}; Surface(22) = {22};
Line Loop(23) = {73, -(21+100), -72, 21}; Surface(23) = {23};
Line Loop(24) = {74, -(22+100), -73, 22}; Surface(24) = {24};
Line Loop(25) = {83, -(31+100), -82, 31}; Surface(25) = {25};
Line Loop(26) = {82, -(30+100), -81, 30}; Surface(26) = {26};

Transfinite Line {11+100,20+100,21+100,22+100,31+100} = N_x_flow;
Transfinite Line {10+100,30+100} = N_x_flow/2;
Transfinite Line {61,62,63,71,72,73,74,81,82,83} = N_z_flow Using Bump R_z_flow;

//Physical Tags
If (Which_Mesh_Part==1)
    //Physical Surface("fluid_interface") = {10, 11, 12, 13, 14, 15, 16, 20, 21, 22, 23, 24, 25, 26};
    Physical Surface("fluid_bottom_interface") = {13, 14, 15, 12, 11, 16, 10};
    Physical Surface("fluid_pin1") = {20, 21};
    Physical Surface("fluid_pin2") = {24, 23, 22};
    Physical Surface("fluid_pin3") = {25, 26};

ElseIf (Which_Mesh_Part==2)
    //Physical Surface("solid_interface") = {13, 14, 15, 12, 11, 16, 10, 20, 21, 24, 23, 22, 25, 26};
    Physical Surface("solid_bottom_interface") = {13, 14, 15, 12, 11, 16, 10};
    Physical Surface("solid_pin1") = {20, 21};
    Physical Surface("solid_pin2") = {24, 23, 22};
    Physical Surface("solid_pin3") = {25, 26};

EndIf
// ------------------------------------------------------------------------- //
If (Which_Mesh_Part == 0 || Which_Mesh_Part == 1) // fluid only
    // Fluid only description
    // upper additional structured mesh points
    Point(40+100) = {length/4 + Tan(30*rad2deg)*width/2, width, upper_h, gs}; // first half, large y
    Point(41+100) = {length/4 - Tan(30*rad2deg)*width/2, 0, upper_h, gs}; // first half, small y
    Point(42+100) = {length*3/4 - Tan(30*rad2deg)*width/2, width, upper_h, gs}; // second half, large y
    Point(43+100) = {length*3/4 + Tan(30*rad2deg)*width/2, 0, upper_h, gs}; // second half, small y
    Point(44+100) = {0, 0, upper_h, gs}; // corner point inlet
    Point(45+100) = {length, 0, upper_h, gs}; // corner point outlet

    // y+ setting helper points
    If(1==0)
        Point(10000) = {0, width-r_pin_upper, upper_h, gs}; // lower pin1 on inlet
        Point(10001) = {0, width-r_pin_upper-2.0e-5, upper_h, gs}; // lower pin1 on inlet
        Line(10000) = {10000,10001};
    EndIf
    // upper additional structured mesh lines
    // outer boundary
    Line(40+100) = {11+100, 44+100};
    Line(41+100) = {44+100, 41+100};
    Line(42+100) = {41+100, 21+100};
    Line(43+100) = {43+100, 24+100};
    Line(44+100) = {43+100, 45+100};
    Line(45+100) = {45+100, 31+100};
    Line(46+100) = {33+100, 42+100};
    Line(47+100) = {42+100, 40+100};
    Line(48+100) = {40+100, 13+100};
    // inner lines
    Line(49+100) = {41+100, 12+100};
    Line(50+100) = {41+100, 40+100};
    Line(51+100) = {22+100, 40+100};
    Line(52+100) = {23+100, 42+100};
    Line(53+100) = {42+100, 43+100};
    Line(54+100) = {43+100, 32+100};

    // line loops and surfaces on lower domain interface
    Line Loop(10+100) = {40+100, 41+100, 49+100, -(10+100)}; Plane Surface(10+100) = {10+100};
    Line Loop(11+100) = {-(49+100), 50+100, 48+100, -(11+100)}; Plane Surface(11+100) = {11+100};
    Line Loop(12+100) = {42+100, 20+100, 51+100, -(50+100)}; Plane Surface(12+100) = {12+100};
    Line Loop(13+100) = {-(51+100), 21+100, 52+100, 47+100}; Plane Surface(13+100) = {13+100};
    Line Loop(14+100) = {53+100, 43+100, -(22+100), 52+100}; Plane Surface(14+100) = {14+100};
    Line Loop(15+100) = {53+100, 54+100, 31+100, 46+100}; Plane Surface(15+100) = {15+100};
    Line Loop(16+100) = {44+100, 45+100, 30+100, -(54+100)}; Plane Surface(16+100) = {16+100};

    // No progression in flow direction on the pin surface
    Transfinite Line {50+100,47+100,53+100} = N_x_flow;
    Transfinite Line {41+100,44+100} = N_x_flow/2;
    // Progression normal to the pin surface
    Transfinite Line {40+100, -(49+100), -(48+100), -(42+100), 51+100, 52+100, -(43+100), 46+100, -(54+100), -(45+100)} = N_y_flow Using Progression R_y_flow;

    // Side Faces on the outside of the domain
    // additional lines in z-direction on fluid boundary
    Line(160) = {44, 144};
    Line(161) = {41, 141};
    Line(162) = {43, 143};
    Line(163) = {45, 145};
    Line(164) = {42, 142};
    Line(165) = {40, 140};
    Transfinite Line {160,161,162,163,164,165} = N_z_flow Using Bump 0.1;

    // Side-faces
    Line Loop(117) = {61, 140, -160, -40}; Surface(117) = {117};
    Line Loop(118) = {160, 141, -161, -41}; Surface(118) = {118};
    Line Loop(119) = {161, 142, -71, -42}; Surface(119) = {119};
    Line Loop(120) = {74, -143, -162, 43}; Surface(120) = {120};
    Line Loop(121) = {162, 144, -163, -44}; Surface(121) = {121};
    Line Loop(122) = {163, 145, -81, -45}; Surface(122) = {122};
    Line Loop(123) = {83, 146, -164, -46}; Surface(123) = {123};
    Line Loop(124) = {164, 147, -165, -47}; Surface(124) = {124};
    Line Loop(125) = {165, 148, -63, -48}; Surface(125) = {125};

    // Internal fluid faces
    Line Loop(126) = {62, -149, -161, 49}; Surface(126) = {126};
    Line Loop(127) = {165, -150, -161, 50}; Surface(127) = {127};
    Line Loop(128) = {165, -151, -72, 51}; Surface(128) = {128};
    Line Loop(129) = {164, -152, -73, 52}; Surface(129) = {129};
    Line Loop(130) = {164, 153, -162, -53}; Surface(130) = {130};
    Line Loop(131) = {82, -154, -162, 54}; Surface(131) = {131};

    // Definition
    Surface Loop(1) = {110, 117, 20, 10, 118, 126}; Volume(1) = {1};
    Surface Loop(2) = {111, 125, 21, 11, 126, 127}; Volume(2) = {2};
    Surface Loop(3) = {112, 119, 22, 12, 127, 128}; Volume(3) = {3};
    Surface Loop(4) = {113, 23, 13, 124, 128, 129}; Volume(4) = {4};
    Surface Loop(5) = {114, 120, 24, 14, 129, 130}; Volume(5) = {5};
    Surface Loop(6) = {115, 25, 123, 15, 130, 131}; Volume(6) = {6};
    Surface Loop(7) = {116, 121, 122, 26, 16, 131}; Volume(7) = {7};

    //Physical Tags
    Physical Volume("fluid_volume")     = {1, 2, 3, 4, 5, 6, 7};
    Physical Surface("fluid_top")       = {110, 111, 112, 113, 114, 115, 116};
    Physical Surface("fluid_inlet")     = {117};
    Physical Surface("fluid_outlet")    = {122};
    Physical Surface("fluid_sym_sides") = {119, 118, 120, 121, 123, 124, 125};

EndIf

// ------------------------------------------------------------------------- //
// Solid only description
If (Which_Mesh_Part == 0 || Which_Mesh_Part == 2)
    
    If (1==1) // Pin 1 Solid
        // Solid inner pin 1 and bottom300-er range
        // pin 1 lower 
        Point(301) = {InnerRadiusFactor*0, width-r_pin_lower*InnerRadiusFactor, 0, gs}; // lower pin1 on inlet
        Point(302) = {InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, width-Cos(30*rad2deg)*r_pin_lower*InnerRadiusFactor, 0, gs}; // lower pin1 in between
        Point(303) = {InnerRadiusFactor*r_pin_lower, width, 0, gs}; // lower pin1 on sym
        Line(301) = {301,302};
        Line(302) = {302,303};

        // pin 1 upper
        Point(304) = {InnerRadiusFactor*0, width-r_pin_upper*InnerRadiusFactor, upper_h, gs}; // lower pin1 on inlet
        Point(305) = {InnerRadiusFactor*Sin(30*rad2deg)*r_pin_upper, width-Cos(30*rad2deg)*r_pin_upper*InnerRadiusFactor, upper_h, gs}; // lower pin1 in between
        Point(306) = {InnerRadiusFactor*r_pin_upper, width, upper_h, gs}; // lower pin1 on sym
        Line(303) = {304,305};
        Line(304) = {305,306};

        // pin 1 additional lines
        Line(305) = {10, 301};
        Line(306) = {301, 11};
        Line(307) = {10, 303};
        Line(308) = {303, 13};
        Line(309) = {110, 304};
        Line(310) = {304, 111};
        Line(311) = {110, 306};
        Line(312) = {306, 113};
        Line(317) = {305, 112};
        Line(318) = {302, 12};

        // pin1 in z-direction
        Line(313) = {301, 304};
        Line(314) = {302, 305};
        Line(315) = {303, 306};
        Line(316) = {10, 110};

        // structured mesh
        Transfinite Line {304,302,305,309} = N_x_flow;
        Transfinite Line {303,301,307,311} = N_x_flow/2;
        Transfinite Line {306,318,308,310,317,312} = N_y_innerPin Using Progression R_y_innerPin;
        Transfinite Line {313,314,315,316} = N_z_flow Using Bump R_z_flow;

        // pin1 Line loop and surface definition
        Line Loop(27) = {313, 310, -61, -306}; Plane Surface(27) = {27};
        Line Loop(28) = {316, 309, -313, -305}; Plane Surface(28) = {28};
        Line Loop(29) = {315, -311, -316, 307}; Plane Surface(29) = {29};
        Line Loop(30) = {315, -304, -314, 302}; Surface(30) = {30};
        Line Loop(31) = {314, -303, -313, 301}; Surface(31) = {31};
        Line Loop(32) = {308, -11, -318, 302}; Plane Surface(32) = {32};
        Line Loop(33) = {318, -10, -306, 301}; Plane Surface(33) = {33};
        Line Loop(34) = {307, -302, -301, -305}; Plane Surface(34) = {34};
        Line Loop(35) = {304, 312, -111, -317}; Plane Surface(35) = {35};
        Line Loop(36) = {317, -110, -310, 303}; Plane Surface(36) = {36};
        Line Loop(37) = {311, -304, -303, -309}; Plane Surface(37) = {37};
        Line Loop(38) = {63, -312, -315, 308}; Plane Surface(38) = {38};
        Line Loop(39) = {317, -62, -318, 314}; Plane Surface(39) = {39};

        Surface Loop(8) = {33, 27, 36, 20, 39, 31}; Volume(8) = {8};
        Surface Loop(9) = {32, 38, 21, 35, 30, 39}; Volume(9) = {9};
        Surface Loop(10) = {37, 29, 28, 30, 31, 34};Volume(10) = {10};

        Physical Volume("solid_volume") = {8,9,10};

        //Physical Tags

    EndIf

    If (1==1) // Pin 2 solid
        // Solid inner half pin 2 and bottome300-er range (copied from interface and solid pin parts)
        // Lower Pin2
        Point(320)  = {0.5*length, 0, 0, gs}; // pin midpoint
        Point(321)  = {0.5*length - InnerRadiusFactor*r_pin_lower, InnerRadiusFactor*0, 0, gs}; // lower small x
        Point(322)  = {length/2 - InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, InnerRadiusFactor*Cos(30*rad2deg)*r_pin_lower, 0, gs}; // small intermediate
        Point(323)  = {length/2 + InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, InnerRadiusFactor*Cos(30*rad2deg)*r_pin_lower, 0, gs}; // large intermediate
        Point(324)  = {0.5*length + InnerRadiusFactor*r_pin_lower, InnerRadiusFactor*0, 0, gs}; // lower large x
        Line(320) = {321,322}; // first segment
        Line(321) = {322,323}; // second segment
        Line(322) = {323,324}; // third segment

        // Upper Pin2
        Point(330)  = {0.5*length, 0, upper_h, gs}; // pin midpoint
        Point(331)  = {0.5*length - InnerRadiusFactor*r_pin_upper, InnerRadiusFactor*0, upper_h, gs}; // lower small x
        Point(332)  = {length/2 - InnerRadiusFactor*Sin(30*rad2deg)*r_pin_upper, InnerRadiusFactor*Cos(30*rad2deg)*r_pin_upper, upper_h, gs}; // small intermediate
        Point(333)  = {length/2 + InnerRadiusFactor*Sin(30*rad2deg)*r_pin_upper, InnerRadiusFactor*Cos(30*rad2deg)*r_pin_upper, upper_h, gs}; // large intermediate
        Point(334)  = {0.5*length + InnerRadiusFactor*r_pin_upper, 0, upper_h, gs}; // lower large x
        Line(330) = {331,332}; // first segment
        Line(331) = {332,333}; // second segment
        Line(332) = {333,334}; // third segment

        // pin 2 additional connecting lines
        Line(333) = {21, 321}; // lower
        Line(334) = {22, 322};
        Line(335) = {23, 323};
        Line(336) = {24, 324};
        Line(337) = {324, 321};

        Line(338) = {121, 331}; // upper
        Line(339) = {122, 332};
        Line(340) = {123, 333};
        Line(341) = {124, 334};
        Line(342) = {334, 331};

        //pin 2 connection z-direction (bottom to top)
        Line(343) = {321, 331};
        Line(344) = {322, 332};
        Line(345) = {323, 333};
        Line(346) = {324, 334};

        // Line loops and surfaces
        Line Loop(132) = {71, 338, -343, -333};  Plane Surface(132) = {132};
        Line Loop(133) = {72, 339, -344, -334};  Plane Surface(134) = {133};
        Line Loop(134) = {330, -344, -320, 343};       Surface(135) = {134};
        Line Loop(136) = {120, 339, -330, -338}; Plane Surface(138) = {136};
        Line Loop(137) = {20, 334, -320, -333};  Plane Surface(139) = {137};
        Line Loop(138) = {344, 331, -345, -321};       Surface(140) = {138};
        Line Loop(139) = {345, -340, -73, 335};  Plane Surface(141) = {139};
        Line Loop(140) = {334, 321, -335, -21};  Plane Surface(142) = {140};
        Line Loop(141) = {339, 331, -340, -121}; Plane Surface(143) = {141};
        Line Loop(142) = {345, 332, -346, -322};       Surface(144) = {142};
        Line Loop(143) = {346, -341, -74, 336};  Plane Surface(145) = {143};
        Line Loop(144) = {332, -341, -122, 340}; Plane Surface(146) = {144};
        Line Loop(145) = {322, -336, -22, 335};  Plane Surface(147) = {145};
        Line Loop(146) = {337, 320, 321, 322};   Plane Surface(148) = {146};
        Line Loop(147) = {342, 330, 331, 332};   Plane Surface(149) = {147};
        Line Loop(148) = {343, -342, -346, 337}; Plane Surface(150) = {148};

        // structured parts
        Transfinite Line {-338, -339, -340, -341, -333, -334, -335, -336} = N_y_innerPin Using Progression R_y_innerPin; // lines pointing into circle midpoint
        Transfinite Line {330, 331, 332, 320, 321, 322, 337, 342} = N_x_flow; // circle arcs
        Transfinite Line {343, 344, 345, 346} = N_z_flow Using Bump R_z_flow; // lines in z-direction

        Surface Loop(11) = {138, 22, 132, 134, 139, 135}; Volume(11) = {11};
        Surface Loop(12) = {143, 23, 140, 141, 134, 142}; Volume(12) = {12};
        Surface Loop(13) = {146, 145, 24, 141, 147, 144}; Volume(13) = {13};
        Surface Loop(14) = {149, 150, 144, 135, 140, 148}; Volume(14) = {14};

        Physical Volume("solid_volume") += {11,12,13,14};
    EndIf

    If (1==1) // Pin 3 solid
        // pin 3 structured
        // lower Pin3
        Point(341) = {length, width-InnerRadiusFactor*r_pin_lower, 0, gs}; // on outlet
        Point(342) = {length-InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, width-InnerRadiusFactor*Cos(30*rad2deg)*r_pin_lower,0, gs};
        Point(343) = {length - InnerRadiusFactor*r_pin_lower, width, 0, gs}; // on sym

        // upper Pin3
        Point(344) = {length, width-InnerRadiusFactor*r_pin_upper, upper_h, gs}; // on outlet
        Point(345) = {length-InnerRadiusFactor*Sin(30*rad2deg)*r_pin_upper, width-InnerRadiusFactor*Cos(30*rad2deg)*r_pin_upper,upper_h, gs};
        Point(346) = {length - InnerRadiusFactor*r_pin_upper, width, upper_h, gs}; // on sym

        // lines: inner pin to interface
        Line(347) = {344, 131};
        Line(348) = {345, 132};
        Line(349) = {346, 133};
        Line(350) = {341, 31};
        Line(351) = {342, 32};
        Line(352) = {343, 33};

        Line(353) = {344, 345}; // inner (not circle)
        Line(354) = {345, 346};
        Line(355) = {341, 342};
        Line(356) = {342, 343};

        Line(357) = {341, 344}; // in positive z-direction
        Line(358) = {342, 345};
        Line(359) = {343, 346};
        Line(360) = {30, 130};

        Line(361) = {341, 30}; //core part of the pin lines
        Line(362) = {30, 343};
        Line(363) = {344, 130};
        Line(364) = {130, 346};

        // structured parts
        Transfinite Line {347,348,349,350,351,352} = N_y_innerPin Using Progression R_y_innerPin; // lines pointing into circle midpoint
        Transfinite Line {353, 355, 362, 364} = N_x_flow/2; // circle arcs
        Transfinite Line {354, 356, 361, 363} = N_x_flow; // circle arcs
        Transfinite Line {357, 358, 359, 360} = N_z_flow Using Bump R_z_flow; // lines in z-direction

        //lineloops and surfaces
        Line Loop(149) = {347, 130, -348, -353}; Plane Surface(151) = {149};
        Line Loop(150) = {354, 349, -131, -348}; Plane Surface(152) = {150};
        Line Loop(151) = {364, -354, -353, 363}; Plane Surface(153) = {151};
        Line Loop(152) = {350, 30, -351, -355};  Plane Surface(154) = {152};
        Line Loop(153) = {356, 352, -31, -351};  Plane Surface(155) = {153};
        Line Loop(154) = {362, -356, -355, 361}; Plane Surface(156) = {154};
        Line Loop(155) = {357, 353, -358, -355}; Plane Surface(157) = {155};
        Line Loop(156) = {358, 354, -359, -356}; Plane Surface(158) = {156};
        Line Loop(157) = {81, -347, -357, 350};  Plane Surface(159) = {157};
        Line Loop(158) = {357, 363, -360, -361}; Plane Surface(160) = {158};
        Line Loop(159) = {360, 364, -359, -362}; Plane Surface(161) = {159};
        Line Loop(160) = {359, 349, -83, -352};  Plane Surface(162) = {160};
        Line Loop(446) = {348, -82, -351, 358}; Plane Surface(446) = {446};

        // surface loops
        Surface Loop(15) = {151, 159, 26, 154, 446, 157};  Volume(15) = {15};
        Surface Loop(16) = {152, 162, 25, 155, 158, 446};  Volume(16) = {16};
        Surface Loop(17) = {153, 161, 160, 156, 157, 158}; Volume(17) = {17};

        // physical suf/vol
        Physical Volume("solid_volume") += {15,16,17};
    EndIf

    If (1==1) // Bottom block solid
        // solid bottom (heater) part 400er range
        // Lower Pin1
        Point(410) = {0, width, -lower_h, gs}; // lower pin1 midpoint
        Point(411) = {0, width-r_pin_lower, -lower_h, gs}; // lower pin1 on inlet
        Point(412) = {Sin(30*rad2deg)*r_pin_lower, width-Cos(30*rad2deg)*r_pin_lower, -lower_h, gs}; // lower pin1 in between
        Point(413) = {r_pin_lower, width, -lower_h, gs}; // lower pin1 on sym
        Circle(410) = {411,410,412}; // lower pin1 smaller first part
        Circle(411) = {412,410,413}; // lower pin1 larger second part

        // Lower Pin2
        Point(420) = {0.5*length, 0, -lower_h, gs}; // pin midpoint
        Point(421) = {0.5*length - r_pin_lower, 0, -lower_h, gs}; // lower small x
        Point(422) = {length/2 - Sin(30*rad2deg)*r_pin_lower, Cos(30*rad2deg)*r_pin_lower, -lower_h, gs}; // small intermediate
        Point(423) = {length/2 + Sin(30*rad2deg)*r_pin_lower, Cos(30*rad2deg)*r_pin_lower, -lower_h, gs}; // large intermediate
        Point(424) = {0.5*length + r_pin_lower, 0, -lower_h, gs}; // lower large x
        Circle(420) = {421,420,422}; // first segment
        Circle(421) = {422,420,423}; // second segment
        Circle(422) = {423,420,424}; // third segment

        // lower Pin3
        Point(430) = {length, width, -lower_h, gs}; // midpoint
        Point(431) = {length, width-r_pin_lower, -lower_h, gs}; // on outlet
        Point(432) = {length-Sin(30*rad2deg)*r_pin_lower, width-Cos(30*rad2deg)*r_pin_lower,-lower_h, gs};
        Point(433) = {length - r_pin_lower, width, -lower_h, gs}; // on sym
        Circle(430) = {431,430,432}; // first segment
        Circle(431) = {432,430,433}; // second segment

        // lower additional structured mesh points
        Point(440) = {length/4 + Tan(30*rad2deg)*width/2, width, -lower_h, gs}; // first half, large y
        Point(441) = {length/4 - Tan(30*rad2deg)*width/2, 0, -lower_h, gs}; // first half, small y
        Point(442) = {length*3/4 - Tan(30*rad2deg)*width/2, width, -lower_h, gs}; // second half, large y
        Point(443) = {length*3/4 + Tan(30*rad2deg)*width/2, 0, -lower_h, gs}; // second half, small y
        Point(444) = {0, 0, -lower_h, gs}; // corner point inlet
        Point(445) = {length, 0, -lower_h, gs}; // corner point outlet

        // lower additional structured mesh lines
        // outer boundary
        Line(440) = {411, 444};
        Line(441) = {444, 441};
        Line(442) = {441, 421};
        Line(443) = {443, 424};
        Line(444) = {443, 445};
        Line(445) = {445, 431};
        Line(446) = {433, 442};
        Line(447) = {442, 440};
        Line(448) = {440, 413};
        // inner lines
        Line(449) = {441, 412};
        Line(450) = {441, 440};
        Line(451) = {422, 440};
        Line(452) = {423, 442};
        Line(453) = {442, 443};
        Line(454) = {443, 432};

        // line loops and surfaces on lower domain interface
        Line Loop(410) = {440, 441, 449, -410}; Plane Surface(410) = {410};
        Line Loop(411) = {-449, 450, 448, -411}; Plane Surface(411) = {411};
        Line Loop(412) = {442, 420, 451, -450}; Plane Surface(412) = {412};
        Line Loop(413) = {-451, 421, 452, 447}; Plane Surface(413) = {413};
        Line Loop(414) = {453, 443, -422, 452}; Plane Surface(414) = {414};
        Line Loop(415) = {453, 454, 431, 446}; Plane Surface(415) = {415};
        Line Loop(416) = {444, 445, 430, -454}; Plane Surface(416) = {416};

        // No progression in flow direction on the pin surface
        Transfinite Line {411,450,420,447,421,422,453,431} = N_x_flow;
        Transfinite Line {410,441,430,444} = N_x_flow/2;
        // Progression normal to the pin surface
        Transfinite Line {440, -449, -448, -442, 451, 452, -443, 446, -454, -445} = N_y_flow Using Progression R_y_flow; 

        // pin inner parts are in the 500-er range
        // pin 1 lower 
        Point(501) = {InnerRadiusFactor*0, width-r_pin_lower*InnerRadiusFactor, -lower_h, gs}; // lower pin1 on inlet
        Point(502) = {InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, width-Cos(30*rad2deg)*r_pin_lower*InnerRadiusFactor, -lower_h, gs}; // lower pin1 in between
        Point(503) = {InnerRadiusFactor*r_pin_lower, width, -lower_h, gs}; // lower pin1 on sym
        Line(501) = {501,502};
        Line(502) = {502,503};

        Line(503) = {503, 413};
        Line(504) = {502, 412};
        Line(505) = {501, 411};
        Line(506) = {501, 410};
        Line(507) = {410, 503};
        // No progression in flow direction on the pin surface
        Transfinite Line {502,506} = N_x_flow;
        Transfinite Line {501,507} = N_x_flow/2;
        Transfinite Line {503,504,505} = N_y_innerPin Using Progression R_y_innerPin;

        Line Loop(417) = {502, 503, -411, -504}; Plane Surface(417) = {417};
        Line Loop(418) = {501, 504, -410, -505}; Plane Surface(418) = {418};
        Line Loop(419) = {506, 507, -502, -501}; Plane Surface(419) = {419};

        // Lower Pin2
        Point(520)  = {0.5*length, 0, 0, gs}; // pin midpoint
        Point(521)  = {0.5*length - InnerRadiusFactor*r_pin_lower, InnerRadiusFactor*0, -lower_h, gs}; // lower small x
        Point(522)  = {length/2 - InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, InnerRadiusFactor*Cos(30*rad2deg)*r_pin_lower, -lower_h, gs}; // small intermediate
        Point(523)  = {length/2 + InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, InnerRadiusFactor*Cos(30*rad2deg)*r_pin_lower, -lower_h, gs}; // large intermediate
        Point(524)  = {0.5*length + InnerRadiusFactor*r_pin_lower, InnerRadiusFactor*0, -lower_h, gs}; // lower large x
        Line(520) = {521,522}; // first segment
        Line(521) = {522,523}; // second segment
        Line(522) = {523,524}; // third segment

        Line(523) = {521, 421};
        Line(524) = {522, 422};
        Line(525) = {523, 423};
        Line(526) = {524, 424};
        Line(527) = {524, 521};

        // structured parts
        Transfinite Line {523,524,525,526} = N_y_innerPin Using Progression R_y_innerPin; // lines pointing into circle midpoint
        Transfinite Line {520,521,522,527} = N_x_flow; // circle arcs

        Line Loop(420) = {523, 420, -524, -520}; Plane Surface(420) = {420};
        Line Loop(421) = {521, 525, -421, -524}; Plane Surface(421) = {421};
        Line Loop(422) = {522, 526, -422, -525}; Plane Surface(422) = {422};
        Line Loop(423) = {527, 520, 521, 522};   Plane Surface(423) = {423};

        // lower Pin3
        Point(541) = {length, width-InnerRadiusFactor*r_pin_lower, -lower_h, gs}; // on outlet
        Point(542) = {length-InnerRadiusFactor*Sin(30*rad2deg)*r_pin_lower, width-InnerRadiusFactor*Cos(30*rad2deg)*r_pin_lower,-lower_h, gs};
        Point(543) = {length - InnerRadiusFactor*r_pin_lower, width, -lower_h, gs}; // on sym

        Line(528) = {541, 542};
        Line(529) = {542, 543};
        Line(530) = {543, 430};
        Line(531) = {430, 541};
        Line(532) = {541, 431};
        Line(533) = {542, 432};
        Line(534) = {543, 433};

        // structured parts
        Transfinite Line {532,533,534} = N_y_innerPin Using Progression R_y_innerPin; // lines pointing into circle midpoint
        Transfinite Line {528,530} = N_x_flow/2; // circle arcs
        Transfinite Line {529,531} = N_x_flow; // circle arcs

        Line Loop(424) = {532, 430, -533, -528}; Plane Surface(424) = {424};
        Line Loop(425) = {533, 431, -534, -529}; Plane Surface(425) = {425};
        Line Loop(426) = {531, 528, 529, 530};   Plane Surface(426) = {426};

        // connecting lines in z-direction
        Line(535) = {10, 410};
        Line(536) = {301, 501};
        Line(537) = {11, 411};
        Line(538) = {44, 444};
        Line(539) = {302, 502};
        Line(540) = {12, 412};
        Line(541) = {303, 503};
        Line(542) = {41, 441};
        Line(543) = {13, 413};
        Line(544) = {21, 421};
        Line(545) = {40, 440};
        Line(546) = {321, 521};
        Line(547) = {22, 422};
        Line(548) = {322, 522};
        Line(549) = {20, 20};
        Line(550) = {20, 20};
        Line(551) = {324, 524};
        Line(552) = {323, 523};
        Line(553) = {23, 423};
        Line(554) = {42, 442};
        Line(555) = {24, 424};
        Line(556) = {43, 443};
        Line(557) = {33, 433};
        Line(558) = {343, 543};
        Line(559) = {32, 432};
        Line(560) = {342, 542};
        Line(561) = {30, 430};
        Line(562) = {341, 541};
        Line(563) = {31, 431};
        Line(564) = {45, 445};
        Transfinite Line {535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564} = N_z_solid Using Progression R_z_solid;

        // lineloops and surfaces of solid inner
        Line Loop(457) = {440, -538, -40, 537};           Plane Surface(457) = {457};
        Line Loop(458) = {538, 441, -542, -41};           Plane Surface(458) = {458};
        Line Loop(459) = {449, -540, -49, 542};           Plane Surface(459) = {459};
        Line Loop(460) = {537, 410, -540, -10};           Surface(460) = {460};
        Surface Loop(28) = {410, 457, 458, 459, 460, 10}; Volume(28) = {28};
        Line Loop(461) = {542, 450, -545, -50};           Plane Surface(461) = {461};
        Line Loop(462) = {448, -543, -48, 545};           Plane Surface(462) = {462};
        Line Loop(463) = {411, -543, -11, 540};           Surface(463) = {463};
        Surface Loop(29) = {411, 462, 463, 11, 461, 459}; Volume(29) = {29};
        Line Loop(464) = {542, 442, -544, -42};           Plane Surface(464) = {464};
        Line Loop(465) = {451, -545, -51, 547};           Plane Surface(465) = {465};
        Line Loop(466) = {547, -420, -544, 20};           Surface(466) = {466};
        Surface Loop(30) = {412, 464, 466, 465, 12, 461}; Volume(30) = {30};
        Line Loop(467) = {447, -545, -47, 554};           Plane Surface(467) = {467};
        Line Loop(468) = {452, -554, -52, 553};           Plane Surface(468) = {468};
        Line Loop(469) = {421, -553, -21, 547};           Surface(469) = {469};
        Surface Loop(31) = {413, 467, 13, 468, 469, 465}; Volume(31) = {31};
        Line Loop(470) = {453, -556, -53, 554};           Plane Surface(470) = {470};
        Line Loop(471) = {443, -555, -43, 556};           Plane Surface(471) = {471};
        Line Loop(472) = {422, -555, -22, 553};           Surface(472) = {472};
        Line Loop(473) = {454, -559, -54, 556};           Plane Surface(473) = {473};
        Surface Loop(32) = {414, 471, 472, 14, 470, 468}; Volume(32) = {32};
        Line Loop(474) = {446, -554, -46, 557};           Plane Surface(474) = {474};
        Line Loop(475) = {431, -557, -31, 559};           Surface(475) = {475};
        Surface Loop(33) = {415, 474, 15, 475, 473, 470}; Volume(33) = {33};
        Line Loop(476) = {444, -564, -44, 556};           Plane Surface(476) = {476};
        Line Loop(477) = {564, 445, -563, -45};           Plane Surface(477) = {477};
        Line Loop(478) = {430, -559, -30, 563};           Surface(478) = {478};
        Surface Loop(34) = {416, 476, 477, 478, 16, 473}; Volume(34) = {34};

        Physical Volume("solid_volume") += {28,29,30,31,32,33,34};

        If (1==1) // pin 1 inside lower block 
            // pin 1
            Line Loop(427) = {535, -506, -536, -305}; Plane Surface(427) = {427};
            Line Loop(428) = {507, -541, -307, 535};  Plane Surface(428) = {428};
            Line Loop(429) = {536, 505, -537, -306};  Plane Surface(429) = {429};
            Line Loop(430) = {503, -543, -308, 541};  Plane Surface(430) = {430};
            Line Loop(431) = {502, -541, -302, 539};  Plane Surface(431) = {431};
            Line Loop(432) = {501, -539, -301, 536};  Plane Surface(432) = {432};
            Line Loop(433) = {543, -411, -540, 11};         Surface(433) = {433};
            Line Loop(434) = {410, -540, -10, 537};         Surface(434) = {434};
            Line Loop(447) = {539, 504, -540, -318};  Plane Surface(447) = {447};

            Surface Loop(18) = {34, 427, 428, 419, 431, 432}; Volume(18) = {18};
            Surface Loop(19) = {33, 418, 429, 434, 447, 432}; Volume(19) = {19};
            Surface Loop(20) = {417, 430, 433, 447, 32, 431}; Volume(20) = {20};

            Physical Volume("solid_volume") += {18,19,20};
        EndIf

        If (1==1) // pin 2 inside lower block 
            // pin 2
            Line Loop(435) = {544, -523, -546, -333}; Plane Surface(435) = {435};
            Line Loop(436) = {546, -527, -551, 337};  Plane Surface(436) = {436};
            Line Loop(437) = {551, 526, -555, 336};   Plane Surface(437) = {437};
            Line Loop(438) = {520, -548, -320, 546};  Plane Surface(438) = {438};
            Line Loop(439) = {548, 521, -552, -321};  Plane Surface(439) = {439};
            Line Loop(440) = {552, 522, -551, -322};  Plane Surface(440) = {440};
            Line Loop(441) = {524, -547, 334, 548};   Plane Surface(441) = {441};
            Line Loop(442) = {553, -525, -552, -335}; Plane Surface(442) = {442};
            Line Loop(443) = {420, -547, -20, 544};         Surface(443) = {443};
            Line Loop(444) = {547, 421, -553, -21};         Surface(444) = {444};
            Line Loop(445) = {422, -555, -22, 553};         Surface(445) = {445};

            // surface loops
            Surface Loop(21) = {139, 435, 443, 420, 438, 441}; Volume(21) = {21};
            Surface Loop(22) = {439, 421, 441, 442, 444, 142}; Volume(22) = {22};
            Surface Loop(23) = {422, 437, 445, 440, 442, 147}; Volume(23) = {23};
            Surface Loop(24) = {423, 436, 438, 439, 440, 148}; Volume(24) = {24};

            Physical Volume("solid_volume") += {21,22,23,24};
        EndIf

        If (1==1) // pin 3 inside lower block
            Line Loop(448) = {531, -562, 361, 561};  Plane Surface(448) = {448};
            Line Loop(449) = {562, 528, -560, -355}; Plane Surface(449) = {449};
            Line Loop(450) = {529, -558, -356, 560}; Plane Surface(450) = {450};
            Line Loop(451) = {558, 530, -561, 362};  Plane Surface(451) = {451};
            Line Loop(452) = {532, -563, -350, 562}; Plane Surface(452) = {452};
            Line Loop(453) = {533, -559, -351, 560}; Plane Surface(453) = {453};
            Line Loop(454) = {534, -557, -352, 558}; Plane Surface(454) = {454};
            Line Loop(455) = {563, 430, -559, -30};        Surface(455) = {455};
            Line Loop(456) = {431, -557, -31, 559};        Surface(456) = {456};

            Surface Loop(25) = {426, 448, 451, 450, 449, 156}; Volume(25) = {25};
            Surface Loop(26) = {424, 452, 455, 453, 154, 449}; Volume(26) = {26};
            Surface Loop(27) = {425, 454, 456, 450, 453, 155}; Volume(27) = {27};

            Physical Volume("solid_volume") += {25,26,27};
        EndIf

    EndIf // Bottom Block solid

    Physical Surface("solid_sym_sides")     = {451, 454, 474, 467, 462, 430, 428, 458, 464, 435, 436, 437, 471, 476, 150, 145, 132, 38, 29, 161, 162};
    Physical Surface("solid_pins_top")      = {37, 36, 35, 138, 143, 149, 146, 152, 153, 151};
    //Physical Surface("solid_pin1_top")      = {37, 36, 35};
    //Physical Surface("solid_pin2_top")      = {138, 143, 149, 146};
    //Physical Surface("solid_pin3_top")      = {152, 153, 151};
    Physical Surface("solid_bottom_heater") = {416, 424, 425, 415, 414, 426, 413, 421, 422, 423, 420, 412, 411, 410, 418, 417, 419};
    Physical Surface("solid_pin3_outlet")   = {159, 160};
    Physical Surface("solid_pin1_inlet")    = {28, 27};
    Physical Surface("solid_block_inlet")   = {427, 429, 457};
    Physical Surface("solid_block_outlet")  = {477, 452, 448};

EndIf

// ------------------------------------------------------------------------- //
// Meshing
Coherence;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

If (Do_Meshing == 1)
    Mesh 1; Mesh 2; Mesh 3;
EndIf

// ------------------------------------------------------------------------- //
// Write .su2 meshfile
If (Write_mesh == 1)

    Mesh.Format = 42; // .su2 mesh format, 
    If (Which_Mesh_Part == 1)
        Save "fluid.su2";
    ElseIf (Which_Mesh_Part == 2)
        Save "solid.su2";
    Else
        Printf("Unvalid Which_Mesh_Part variable for output writing.");
        Abort;
    EndIf

EndIf

// Write .cgns meshfile
//If (Write_mesh == 1)
//
//    Mesh.Format = 32; // .cgns mesh format, 
//    If (Which_Mesh_Part == 1)
//        Save "fluid.cgns";
//    ElseIf (Which_Mesh_Part == 2)
//        Save "solid.cgns";
//    ElseIf (Which_Mesh_Part == 0)
//        Save "fluid_and_solid.cgns";
//    Else
//        Printf("Unvalid Which_Mesh_Part variable for output writing.");
//        Abort;
//    EndIf
//
//EndIf

