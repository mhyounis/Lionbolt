// Gmsh project created on Sun May 21 18:58:20 2023
SetFactory("OpenCASCADE");

// DESCRIPTION
  // This is a mesh for a single-material problem wherein a box is being irradiated by an external beam which is either a planar or spherical source, depending on the user-input parameters below. This is good for water tank calculations, after which the model is named. However, it can be used for any material as long as the material is defined within wiscobolt. The box is sitting on the origin of the xy-plane. That is to say, the center of mass is always (0,0,Lz/2). The mesh has planes at x = 0 and y = 0 which are useful for visualization.

// USER INPUT
  // Mesh element size (can also be changed in GMSH program)
    lc = 0.1;

  // Scale
    Lx = 1;   // Size of box along x
    Ly = 1;   // Size of box along y
    Lz = 1;   // Size of box along z
    d1 = 0.1; // Fraction of Lz at which the slab starts (facing down)
    d2 = 0.2; // Fraction of Lz at which the slab ends   (facing down)

  // Beam profile
  // Notes: For a spherical source, fsbot > fstop. For a planar source, you can set fstop = fsbot. There is no significance to fsbot < fstop, and my program will not provide physically meaningful results for this case.
  //        It is your responsibility to make sure that the appropriate settings/apex angles are input to the program when this mesh is being used, as it will not be checked automatically based on the mesh.
    fsxtop = 1/3; // Width of beam cutout at z = Lz along x
    fsytop = 1/3; // Width of beam cutout at z = Lz along y
    fsxbot = (2/5)*(1+5/60); // Width of beam cutout at z = 0 along x
    fsybot = (2/5)*(1+5/60); // Width of beam cutout at z = 0 along y

// GEOMETRY
  // (The following can be ignored)
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {Lx, 0, 0, lc};
//+
Point(3) = {0, Ly, 0, lc};
//+
Point(4) = {0, 0, Lz, lc};
//+
Point(5) = {Lx, Ly, 0, lc};
//+
Point(6) = {Lx, 0, Lz, lc};
//+
Point(7) = {0, Ly, Lz, lc};
//+
Point(8) = {Lx, Ly, Lz, lc};
//+
Point(9) = {Lx/2-fsxtop/2, Ly/2-fsytop/2, Lz, lc};
//+
Point(10) = {Lx/2+fsxtop/2, Ly/2-fsytop/2, Lz, lc};
//+
Point(11) = {Lx/2-fsxtop/2, Ly/2+fsytop/2, Lz, lc};
//+
Point(12) = {Lx/2+fsxtop/2, Ly/2+fsytop/2, Lz, lc};
//+
Point(13) = {Lx/2+fsxbot/2, Ly/2+fsybot/2, 0, lc};
//+
Point(14) = {Lx/2+fsxbot/2, Ly/2-fsybot/2, 0, lc};
//+
Point(15) = {Lx/2-fsxbot/2, Ly/2-fsybot/2, 0, lc};
//+
Point(16) = {Lx/2-fsxbot/2, Ly/2+fsybot/2, 0, lc};
//+
Point(17) = {Lx/2, 0, 0, lc};
//+
Point(18) = {Lx/2, 0, Lz, lc};
//+
Point(19) = {Lx/2, Ly, Lz, lc};
//+
Point(20) = {Lx/2, Ly, 0, lc};
//+
Point(21) = {0, Ly/2, 0, lc};
//+
Point(22) = {0, Ly/2, Lz, lc};
//+
Point(23) = {Lx, Ly/2, Lz, lc};
//+
Point(24) = {Lx, Ly/2, 0, lc};
//+
Point(25) = {Lx/2, Ly/2, Lz, lc};
//+
Point(26) = {Lx/2, Ly/2, 0, lc};
//+
Point(27) = {Lx/2, Ly/2-fsytop/2, Lz, lc};
//+
Point(28) = {Lx/2, Ly/2+fsytop/2, Lz, lc};
//+
Point(29) = {Lx/2-fsxtop/2, Ly/2, Lz, lc};
//+
Point(30) = {Lx/2+fsxtop/2, Ly/2, Lz, lc};
//+
Point(31) = {Lx/2+fsxbot/2, Ly/2, 0, lc};
//+
Point(32) = {Lx/2, Ly/2+fsybot/2, 0, lc};
//+
Point(33) = {Lx/2-fsxbot/2, Ly/2, 0, lc};
//+
Point(34) = {Lx/2, Ly/2-fsybot/2, 0, lc};
//+
Line(1) = {11, 29};
//+
Line(2) = {29, 25};
//+
Line(3) = {25, 28};
//+
Line(4) = {28, 11};
//+
Line(5) = {28, 12};
//+
Line(6) = {12, 30};
//+
Line(7) = {30, 25};
//+
Line(8) = {25, 27};
//+
Line(9) = {27, 10};
//+
Line(10) = {10, 30};
//+
Line(11) = {29, 9};
//+
Line(12) = {9, 27};
//+
Line(13) = {27, 18};
//+
Line(14) = {18, 6};
//+
Line(15) = {6, 23};
//+
Line(16) = {23, 30};
//+
Line(17) = {8, 23};
//+
Line(18) = {8, 19};
//+
Line(19) = {19, 28};
//+
Line(20) = {19, 7};
//+
Line(21) = {7, 22};
//+
Line(22) = {22, 29};
//+
Line(23) = {22, 4};
//+
Line(24) = {4, 18};
//+
Line(25) = {7, 3};
//+
Line(26) = {3, 20};
//+
Line(27) = {20, 19};
//+
Line(28) = {20, 5};
//+
Line(29) = {5, 8};
//+
Line(30) = {23, 24};
//+
Line(31) = {24, 5};
//+
Line(32) = {24, 2};
//+
Line(33) = {2, 6};
//+
Line(34) = {18, 17};
//+
Line(35) = {17, 2};
//+
Line(36) = {26, 34};
//+
Line(37) = {34, 14};
//+
Line(38) = {14, 31};
//+
Line(39) = {31, 26};
//+
Line(40) = {26, 32};
//+
Line(41) = {32, 13};
//+
Line(42) = {13, 31};
//+
Line(43) = {32, 16};
//+
Line(44) = {16, 33};
//+
Line(45) = {33, 26};
//+
Line(46) = {33, 15};
//+
Line(47) = {15, 34};
//+
Line(48) = {17, 1};
//+
Line(49) = {1, 21};
//+
Line(50) = {21, 33};
//+
Line(51) = {21, 3};
//+
Line(52) = {22, 21};
//+
Line(53) = {4, 1};
//+
Line(54) = {32, 20};
//+
Line(55) = {34, 17};
//+
Line(56) = {31, 24};
//+
Line(57) = {10, 14};
//+
Line(58) = {12, 13};
//+
Line(59) = {28, 32};
//+
Line(60) = {11, 16};
//+
Line(61) = {29, 33};
//+
Line(62) = {15, 9};
//+
Line(63) = {30, 31};
//+
Line(64) = {27, 34};
//+
Line(65) = {25, 26};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 7, 8, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, -8, -2, 11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {19, 4, 1, -22, -21, -20};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {18, 19, 5, 6, -16, -17};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {15, 16, -10, -9, 13, 14};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {22, 11, 12, 13, -24, -23};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {44, 45, 40, 43};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {41, 42, 39, 40};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {38, 39, 36, 37};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {47, -36, -45, 46};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {51, 26, -54, 43, 44, -50};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {28, -31, -56, -42, -41, 54};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {32, -35, -55, 37, 38, 56};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {55, 48, 49, 50, 46, 47};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {24, 34, 48, -53};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {34, 35, 33, -14};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {33, 15, 30, 32};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {17, 30, 31, 29};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {28, 29, 18, -27};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {27, 20, 25, 26};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {25, -51, -52, -21};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {49, -52, 23, 53};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {64, 37, -57, -9};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {64, -36, -65, 8};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {65, -39, -63, 7};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {38, -63, -10, 57};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {58, 42, -63, -6};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {5, 58, -41, -59};
//+
Plane Surface(30) = {30};
//+
Curve Loop(31) = {59, -40, -65, 3};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {4, 60, -43, -59};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {60, 44, -61, -1};
//+
Plane Surface(33) = {33};
//+
Curve Loop(34) = {61, 45, -65, -2};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {11, -62, -46, -61};
//+
Plane Surface(35) = {35};
//+
Curve Loop(36) = {62, 12, 64, -47};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {13, 34, -55, -64};
//+
Plane Surface(37) = {37};
//+
Curve Loop(38) = {56, -30, 16, 63};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {19, 59, 54, 27};
//+
Plane Surface(39) = {39};
//+
Curve Loop(40) = {22, 61, -50, -52};
//+
Plane Surface(40) = {40};
//+
Surface Loop(1) = {26, 27, 28, 25, 3, 11};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {2, 31, 30, 29, 10, 27};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {1, 33, 32, 9, 34, 31};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {34, 35, 36, 26, 4, 12};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {8, 17, 16, 24, 40, 36, 35, 37};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {38, 28, 25, 37, 15, 19, 18, 7};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {20, 6, 21, 14, 29, 30, 39, 38};
//+
Volume(7) = {7};
//+
Surface Loop(8) = {39, 22, 5, 23, 13, 40, 33, 32};
//+
Volume(8) = {8};
//+
Point(35) = {0, 0, Lz - d1, lc};
//+
Point(36) = {0, 0, Lz - d2, lc};
//+
Point(37) = {Lx, 0, Lz - d1, lc};
//+
Point(38) = {Lx, 0, Lz - d2, lc};
//+
Point(39) = {0, Ly, Lz - d1, lc};
//+
Point(40) = {0, Ly, Lz - d2, lc};
//+
Point(41) = {Lx, Ly, Lz - d1, lc};
//+
Point(42) = {Lx, Ly, Lz - d2, lc};
//+
Line(66) = {35, 36};
//+
Line(67) = {36, 38};
//+
Line(68) = {38, 37};
//+
Line(69) = {37, 35};
//+
Line(70) = {37, 41};
//+
Line(71) = {41, 42};
//+
Line(72) = {42, 38};
//+
Line(73) = {42, 40};
//+
Line(74) = {40, 39};
//+
Line(75) = {39, 41};
//+
Line(76) = {39, 35};
//+
Line(77) = {36, 40};
//+
Curve Loop(41) = {66, 67, 68, 69};
//+
Plane Surface(41) = {41};
//+
Curve Loop(42) = {68, 70, 71, 72};
//+
Plane Surface(42) = {42};
//+
Curve Loop(43) = {73, 74, 75, 71};
//+
Plane Surface(43) = {43};
//+
Curve Loop(44) = {77, 74, 76, 66};
//+
Plane Surface(44) = {44};
//+
Curve Loop(45) = {77, -73, 72, -67};
//+
Plane Surface(45) = {45};
//+
Curve Loop(46) = {70, -75, 76, -69};
//+
Plane Surface(46) = {46};
//+
Surface Loop(9) = {45, 44, 43, 46, 42, 41};
//+
Volume(9) = {9};
//+
BooleanFragments{ Volume{6}; Volume{5}; Volume{1}; Volume{4}; Volume{2}; Volume{3}; Volume{7}; Volume{8}; Delete; }{ Volume{9}; Delete; }
//+
MeshSize {:} = lc;
//+
Field[1] = Box;
//+
Field[1].Thickness = 1;
//+
Field[1].VIn = 0.02;
//+
Field[1].XMax = Lx/2+fsxbot/2;
//+
Field[1].XMin = Lx/2-fsxbot/2;
//+
Field[1].YMax = Ly/2+fsybot/2;
//+
Field[1].YMin = Ly/2-fsybot/2;
//+
Field[1].ZMax = Lz;
//+
Field[1].ZMin = Lz/2;
//+
Background Field = 1;

//+
Physical Volume("mat_1", 1) = {1, 6, 9, 12, 15, 18, 21, 24, 4, 3, 10, 7, 16, 13, 22, 19};
//+
Physical Volume("mat_2", 2) = {5, 2, 11, 8, 17, 14, 23, 20};
