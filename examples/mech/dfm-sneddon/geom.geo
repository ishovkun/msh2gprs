// pi = 3.1415926;
alpha = 30 * Pi / 180;

/* domain_size = 200; // outer */
size_x = 100;
size_y = 200;

lfrac = 5;
thickness = 10;

nf = 32;
h = lfrac / nf;   // fracture element size - width

nx = size_x / h;
ny = size_y / h;
nz = 1;

If (nx % 2 != 0)
    nx = nx + 1;
EndIf
If (ny % 2 != 0)
    ny = ny + 1;
EndIf

left  = 0;
right = size_x;
front = -size_y / 2;
back   = +size_y / 2;

// outer domain boundary
Point(1) = {left,  back,    0, h};  // left top
Point(2) = {left,  front, 0, h};  // left front
Point(3) = {right, front, 0, h};  // right front
Point(4) = {right, back,    0, h};  // right top

// on boundaries
Point(5) = {left, 0, 0, h};  /* left */
Point(6) = {left+lfrac, front, 0, h};  /* front */
Point(7) = {right, 0, 0, h};  /* right */
Point(8) = {left+lfrac, back, 0, h};  /* back */

// frac tip
Point(9) = {left + lfrac, 0, 0, h};

/* outer domain back part */
Line(1) = {1, 5}; // Left
Line(2) = {5, 9}; // frac
Line(3) = {9, 7}; // front right (after tip)
Line(4) = {7, 4}; // right
Line(5) = {4, 8};
Line(6) = {8, 1};
Line(7) = {8, 9};
// outer front part
Line(8) = {5, 2}; // left
Line(9) = {2, 6};
Line(10) = {6, 3};
Line(11) = {3, 7};  /* right */
Line(12) = {6, 9};

// upper
Transfinite Line{1}  = ny/2 + 1;
Transfinite Line{2}  = nf + 1;
Transfinite Line{3}  = (nx - nf) + 1;
Transfinite Line{4}  = ny/2 + 1;
Transfinite Line{5}  = (nx - nf) + 1;
Transfinite Line{6}  = nf + 1;
Transfinite Line{7}  = ny/2 + 1;

// lower
Transfinite Line{8}  = ny/2 + 1;
Transfinite Line{9}  = nf + 1;
Transfinite Line{10}  = (nx - nf) + 1;
Transfinite Line{11}  = ny/2 + 1;
Transfinite Line{12}  = ny/2 + 1;

/* /\* Printf("%g ", (inner_size/2) / hf) ; *\/ */
/* /\* Printf("%g ", (inner_size/2) / hf) ; *\/ */

// LOOPS
// upper
Line Loop(1) = {1, 2, -7, 6};
Line Loop(2) = {7, 3, 4, 5};
// lower
Line Loop(3) = {8, 9, 12, -2};
Line Loop(4) = {-12, 10, 11, -3};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

/* // upper inner */
Transfinite Surface(1) = {1, 5, 9, 8};
Transfinite Surface(2) = {8, 9, 7, 4};
Transfinite Surface(3) = {5, 2, 6, 9};
Transfinite Surface(4) = {9, 6, 3, 7};

Recombine Surface(1);
Recombine Surface(2);
Recombine Surface(3);
Recombine Surface(4);

out1[] = Extrude{0,0,-thickness} { Surface{1}; Layers{1}; Recombine;};
out2[] = Extrude{0,0,-thickness} { Surface{2}; Layers{1}; Recombine;};
out3[] = Extrude{0,0,-thickness} { Surface{3}; Layers{1}; Recombine;};
out4[] = Extrude{0,0,-thickness} { Surface{4}; Layers{1}; Recombine;};

// Labels
Left = 1111111;
Right = 1111112;
Front = 2222221;
Back = 2222222;
Bottom = 3333331;
Top = 3333332;

Frac1 = 1;

Physical Surface(Bottom) = {34, 78, 56, 100};
Physical Surface(Top) = {1, 2, 3, 4};
Physical Surface(Left) = {65, 21};
Physical Surface(Right) = {95, 51};
Physical Surface(Front) = {69, 91};
Physical Surface(Back) = {55, 33};

Physical Surface(Frac1) = {25};

ELASTIC_DRIVER = 9999991;
SDA_DRIVER = 9999992;

Physical Volume(SDA_DRIVER) = {1, 2, 3, 4};

Mesh.Algorithm = 8;
