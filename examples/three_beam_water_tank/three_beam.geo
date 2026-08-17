// Gmsh project created on Mon Aug 03 22:14:27 2026
SetFactory("OpenCASCADE");

Lx = 1.0;
Ly = 1.0;
Lz = 1.0;
fsxtop = 1/3;
fsytop = 1/3;
fsxbot = (2/5)*(1+5/60);
fsybot = (2/5)*(1+5/60);
SSD    = 10/3;

fsxover = (SSD + 2*Lz) * fsxbot / (SSD + Lz);
fsyover = (SSD + 2*Lz) * fsybot / (SSD + Lz);

// Base tank
CUBE1 = newv; Box(CUBE1) = {-Lx/2, -Ly/2, 0, Lx, Ly, Lz};

// Beam 1
p1 = newp; Point(p1) = {-fsxover/2, -fsyover/2, -Lz};
p2 = newp; Point(p2) = { fsxover/2, -fsyover/2, -Lz};
p3 = newp; Point(p3) = { fsxover/2,  fsyover/2, -Lz};
p4 = newp; Point(p4) = {-fsxover/2,  fsyover/2, -Lz};
p5 = newp; Point(p5) = { 	 0,          0, SSD + Lz};

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

// Beam 2 --- Just rotate beam 1. Center the
BEAM2[] = Rotate {{1,0,0}, {0,0,Lz/2}, 4*Pi/3} {Duplicata{Volume{BEAM1};}};


// Beam 3 --- Just rotate beam 2
BEAM3[] = Rotate {{1,0,0}, {0,0,Lz/2}, 8*Pi/3} {Duplicata{Volume{BEAM1};}};

beams[] = {BEAM1, BEAM2, BEAM3};

// 
beams_inside[] = BooleanIntersection{Volume{beams[0], beams[1], beams[2]}; Delete;}{Volume{CUBE1};};
frag[] = BooleanFragments{Volume{CUBE1, beams_inside[]};Delete;}{};

//+

Physical Volume("mat_1", 1) = {15, 6, 12, 19, 5, 11, 10, 17, 13, 9, 21, 18, 7, 16, 2, 4, 14, 8, 20, 3};

//+
Field[1] = Distance;
Field[1].SurfacesList = {96};

Field[2] = Distance;
Field[2].SurfacesList = {85};

Field[3] = Distance;
Field[3].SurfacesList = {58};

Field[4] = Distance;
Field[4].SurfacesList = {42};

Field[5] = Distance;
Field[5].SurfacesList = {43};

Field[6] = Distance;
Field[6].SurfacesList = {45};

Field[7] = Distance;
Field[7].SurfacesList = {47};

Field[8] = Distance;
Field[8].SurfacesList = {38};

Field[9] = Distance;
Field[9].SurfacesList = {30};

// // // // // // // // //

smin = 0.02;
smax = 0.2;
dmin = 0.0;
dmax = 1.0;

Field[10] = Threshold;
Field[10].InField = 1;
Field[10].SizeMin = smin;
Field[10].SizeMax = smax;
Field[10].DistMin = dmin;
Field[10].DistMax = dmax;

Field[11] = Threshold;
Field[11].InField = 2;
Field[11].SizeMin = smin;
Field[11].SizeMax = smax;
Field[11].DistMin = dmin;
Field[11].DistMax = dmax;

Field[12] = Threshold;
Field[12].InField = 3;
Field[12].SizeMin = smin;
Field[12].SizeMax = smax;
Field[12].DistMin = dmin;
Field[12].DistMax = dmax;

Field[13] = Threshold;
Field[13].InField = 4;
Field[13].SizeMin = smin;
Field[13].SizeMax = smax;
Field[13].DistMin = dmin;
Field[13].DistMax = dmax;

Field[14] = Threshold;
Field[14].InField = 5;
Field[14].SizeMin = smin;
Field[14].SizeMax = smax;
Field[14].DistMin = dmin;
Field[14].DistMax = dmax;

Field[15] = Threshold;
Field[15].InField = 6;
Field[15].SizeMin = smin;
Field[15].SizeMax = smax;
Field[15].DistMin = dmin;
Field[15].DistMax = dmax;

Field[16] = Threshold;
Field[16].InField = 7;
Field[16].SizeMin = smin;
Field[16].SizeMax = smax;
Field[16].DistMin = dmin;
Field[16].DistMax = dmax;

Field[17] = Threshold;
Field[17].InField = 8;
Field[17].SizeMin = smin;
Field[17].SizeMax = smax;
Field[17].DistMin = dmin;
Field[17].DistMax = dmax;

Field[18] = Threshold;
Field[18].InField = 9;
Field[18].SizeMin = smin;
Field[18].SizeMax = smax;
Field[18].DistMin = dmin;
Field[18].DistMax = dmax;

Field[19] = Min;
Field[19].FieldsList = {10, 11, 12, 13, 14, 15, 16, 17, 18};
Background Field = 19;
