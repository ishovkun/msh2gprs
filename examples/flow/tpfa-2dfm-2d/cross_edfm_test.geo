// case with 2 fracs

l1 = 160;  // len of frac 1
l2 = 160;  // len of frac 2

a = 80;													// frac 1 angle (deg)
b = 80;													// frac 2 angle (deg)
// Thickness = l / 2;
size_x = 200; // domain size
size_y = 200; // domain size
size_z = 10;

lc = size_x / 10;  /* element size */

//25
x[1] = -size_x / 2;
x[2] = -Sin((90-a)*Pi/180.) * 2*l1;
x[3] = -Sin((90-a)*Pi/180.) * l1/2;
x[4] = -x[3];
x[5] = -x[2];
x[6] = -x[1];

y[1] = -12*l1/2;
y[2] = -Cos((90-a)*Pi/180.) * 2*l1;
y[3] = -Cos((90-a)*Pi/180.) * l1/2;
y[4] = -y[3];
y[5] = -y[2];
y[6] = -y[1];

// frac 1
i=1;
p = newp; Point(p) = {x[3],y[3],-size_z/2, lc}; point[i] = p; i=i+1;
p = newp; Point(p) = {x[4],y[4],-size_z/2, lc}; point[i] = p; i=i+1;

dx = -Sin(a*Pi/180.) * 3*l1;
dy = -Cos(a*Pi/180.) * 3*l1;

// center 3
p = newp; Point(p) = {0,0,-size_z/2, lc}; point[i] = p; i=i+1;

// frac 2
x[3] = Sin((90-b)*Pi/180.) * l2/2;
x[4] = -x[3];
y[3] = -Cos((90-b)*Pi/180.) * l2/2;
y[4] = -y[3];

p = newp; Point(p) = {x[3], y[3], -size_z/2, lc}; point[i] = p; i=i+1;  /* 4 */
p = newp; Point(p) = {x[4], y[4], -size_z/2, lc}; point[i] = p; i=i+1;  /* 5 */

xx[3] = -Sin((90-a)*Pi/180.) * l1/2 - x[3];
xx[4] = -xx[3];
yy[3] = -Cos((90-a)*Pi/180.) * l1/2 - y[3];
yy[4] = -yy[3];

// outer domain
// left front
p = newp; Point(p) = {-size_x/2, -size_y/2, -size_z/2, lc}; point[i] = p; i=i+1;  /* 6 */
// front center
p = newp; Point(p) = {0,         -size_y/2, -size_z/2, lc}; point[i] = p; i=i+1; /* 7 */
// right front
p = newp; Point(p) = {size_x/2,  -size_y/2, -size_z/2, lc}; point[i] = p; i=i+1; /* 8 */
// left center
p = newp; Point(p) = {-size_x/2, 0,         -size_z/2, lc}; point[i] = p; i=i+1;  /* 9 */
// right center
p = newp; Point(p) = {size_x/2,  0,         -size_z/2, lc}; point[i] = p; i=i+1;  /* 10 */
// back left
p = newp; Point(p) = {-size_x/2, size_y/2,  -size_z/2, lc}; point[i] = p; i=i+1;  /* 11 */
// back center
p = newp; Point(p) = {0,         size_y/2,  -size_z/2, lc}; point[i] = p; i=i+1;  /* 12 */
// back right
p = newp; Point(p) = {size_x/2,  size_y/2,  -size_z/2, lc}; point[i] = p; i=i+1;  /* 13 */
//+
Line(1) = {6, 9};
//+
Line(2) = {9, 11};
//+
Line(3) = {11, 12};
//+
Line(4) = {12, 13};
//+
Line(5) = {13, 10};
//+
Line(6) = {10, 8};
//+
Line(7) = {8, 7};
//+
Line(8) = {6, 7};
//+ frac 2 left
Line(9) = {5, 3};
//+ frac 2 right
Line(10) = {3, 4};
/* Transfinite Line{9} = n_frac_vert/2; */
/* Transfinite Line{10} = n_frac_vert/2; */
//+ frac 1 right
Line(11) = {3, 2};
//+ frac 1 left
Line(12) = {1, 3};
/* Transfinite Line{11} = n_frac_vert/2 * l1 / l2; */
/* Transfinite Line{12} = n_frac_vert/2  * l1 / l2; */
//+
Line(13) = {9, 5};
//+
Line(14) = {9, 1};
//+
Line(15) = {2, 10};
//+
Line(16) = {4, 10};
//+
Line(17) = {1, 7};
//+
Line(18) = {7, 4};
//+
Line(19) = {5, 12};
//+
Line(20) = {2, 12};
// +
Line Loop(21) = {17, 18, -10, -12};
//+
Plane Surface(22) = {21};
//+
Line Loop(23) = {10, 16, -15, -11};
//+
Plane Surface(24) = {23};
//+
Line Loop(25) = {9, 11, 20, -19};
//+
Plane Surface(26) = {25};
//+
Line Loop(27) = {14, 12, -9, -13};
//+
Plane Surface(28) = {27};
//+
Line Loop(29) = {2, 3, -19, -13};
//+
Plane Surface(30) = {29};
//+
Line Loop(31) = {1, 14, 17, -8};
//+
Plane Surface(32) = {31};
//+
Line Loop(33) = {18, 16, 6, 7};
//+
Plane Surface(34) = {33};
//+
Line Loop(35) = {15, -5, -4, -20};
//+
Plane Surface(36) = {35};

Recombine Surface(22);
Recombine Surface(24);
Recombine Surface(26);
Recombine Surface(28);
Recombine Surface(30);
Recombine Surface(32);
Recombine Surface(34);
Recombine Surface(36);

vol[] = Extrude{0,0,size_z}{Surface{22,24,26,28,30,32,34,36};Layers{1};Recombine;};

// external boundaries
// Physical Surface(1111111) = {155, 133};
// Physical Surface(1111112) = {185, 203};
// Physical Surface(2222221) = {167, 189};
// Physical Surface(2222222) = {137, 207};
// Physical Surface(3333331) = {34, 32, 36, 30, 24, 22, 28, 26};
// Physical Surface(3333332) = {168, 146, 190, 212, 58, 124, 102, 80};

//  internal boundaries
// Physical Surface(2) = {89, 53};
Physical Surface(1) = {79, 57};

Physical Volume(9999991) = {vol[]};


// Geometry.Tolerance = size_x / lc / 1e3;
Geometry.AutoCoherence = 2;
Mesh.MshFileVersion = 2.2;
// Mesh.MshFileVersion = 4.1;

Mesh 3;  // Generalte 3D mesh
Coherence Mesh;  // Remove duplicate entities
Save "geom.msh";  // Save mesh in MSH format
