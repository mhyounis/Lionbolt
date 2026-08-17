// Gmsh project created on Sun May 21 18:58:20 2023
SetFactory("OpenCASCADE");

// DESCRIPTION
  // This is a mesh for a single-material problem wherein a box is being irradiated by an external beam which is either a planar or spherical source, depending on the user-input parameters below. However, it can be used for any material as long as the material is defined within wiscobolt. The box is sitting on the origin of the xy-plane. That is to say, the center of mass is always (0,0,Lz/2). The mesh has planes at x = 0 and y = 0 which are useful for visualization.

// USER INPUT
  // Mesh element size (can also be changed in GMSH program)
    lc = 0.25;

  // Scale
    Lx = 1; // Size of box along x
    Ly = 1; // Size of box along y
    Lz = 1; // Size of box along z

  // Beam profile
  // Notes: For a spherical source, fsbot > fstop. For a planar source, you can set fstop = fsbot. There is no significance to fsbot < fstop, and my program will not provide physically meaningful results for this case.
  //        It is your responsibility to make sure that the appropriate settings/apex angles are input to the program when this mesh is being used, as it will not be checked automatically based on the mesh.
    S      = 10/3; // Source-to-surface distance
    fsxtop = 1/3; // Width of beam cutout at z = Lz along x
    fsytop = 1/3; // Width of beam cutout at z = Lz along y
    fsxbot = (S + Lz) * fsxtop / S; // Width of beam cutout at z = 0 along x
    fsybot = (S + Lz) * fsytop / S; // Width of beam cutout at z = 0 along y
    fsxover = (S + 2 * Lz) * fsxbot / (S + Lz); // used for beam definition
    fsyover = (S + 2 * Lz) * fsybot / (S + Lz); // user for beam definition

// GEOMETRY

CUBE1 = newv; Box(CUBE1) = {-Lx/2, -Ly/2, 0, Lx, Ly, Lz};

// Beam
p1 = newp; Point(p1) = {-fsxover/2, -fsyover/2, -Lz};
p2 = newp; Point(p2) = { fsxover/2, -fsyover/2, -Lz};
p3 = newp; Point(p3) = { fsxover/2,  fsyover/2, -Lz};
p4 = newp; Point(p4) = {-fsxover/2,  fsyover/2, -Lz};
p5 = newp; Point(p5) = { 	 0,          0, S + Lz};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};
l5 = newl; Line(l5) = {p5,p1};
l6 = newl; Line(l6) = {p5,p2};
l7 = newl; Line(l7) = {p5,p3};
l8 = newl; Line(l8) = {p5,p4};

ll1 = newll; Curve Loop(ll1) = {l1,l2,l3,l4};
s1  = news;  Plane Surface(s1) = {ll1};
ll2 = newll; Curve Loop(ll2) = {l5,l6,l1};
s2  = news;  Plane Surface(s2) = {ll2};
ll3 = newll; Curve Loop(ll3) = {l6,l7,l2};
s3  = news;  Plane Surface(s3) = {ll3};
ll4 = newll; Curve Loop(ll4) = {l7,l8,l3};
s4  = news;  Plane Surface(s4) = {ll4};
ll5 = newll; Curve Loop(ll5) = {l8,l5,l4};
s5  = news;  Plane Surface(s5) = {ll5};

sl1 = newsl; Surface Loop(sl1) = {s1, s2, s3, s4, s5};

BEAM1 = newv;  Volume(BEAM1) = {sl1};

// 
beams_inside[] = BooleanIntersection{Volume{BEAM1}; Delete;}{Volume{CUBE1};};
frag[] = BooleanFragments{Volume{CUBE1, beams_inside[]};Delete;}{};


Field[1] = Box;
//+
Field[1].Thickness = 1;
//+
Field[1].VIn = 0.025;
//+
Field[1].XMax =  fsxbot/2;
//+
Field[1].XMin = -fsxbot/2;
//+
Field[1].YMax =  fsybot/2;
//+
Field[1].YMin = -fsybot/2;
//+
Field[1].ZMax = Lz;
//+
Field[1].ZMin = Lz/2;
//+
Field[2] = Box;
//+
Field[2].Thickness = 1;
//+
Field[2].VIn = 0.04;
//+
Field[2].XMax =  fsxbot/2;
//+
Field[2].XMin = -fsxbot/2;
//+
Field[2].YMax =  fsybot/2;
//+
Field[2].YMin = -fsybot/2;
//+
Field[2].ZMax = Lz;
//+
Field[2].ZMin = 0;
//+
Field[3] = Box;
//+
Field[3].Thickness = 1;
//+
Field[3].VIn = 0.0125;
//+
Field[3].XMax =  fsxbot/2;
//+
Field[3].XMin = -fsxbot/2;
//+
Field[3].YMax =  fsybot/2;
//+
Field[3].YMin = -fsybot/2;
//+
Field[3].ZMax = Lz;
//+
Field[3].ZMin = 9 * Lz / 10;
//+
Field[4] = Min;
//+
Field[4].FieldsList = {1, 2, 3};
//+
Background Field = 4;
//+
Physical Volume("mat_1", 1) = {33, 32};
