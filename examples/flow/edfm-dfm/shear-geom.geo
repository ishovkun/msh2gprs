// pi = 3.1415926;
alpha = 30 * Pi / 180;

domain_size = 200; // outer
lfrac = 30;
inner_size = 4*lfrac;
thickness = -10;
nz = 1;

nf = 4;
h = lfrac / nf;  // outer boundary element size
hf = lfrac / nf;   // fracture element size - width

left  = - domain_size / 2;
right = + domain_size / 2;
front = - domain_size / 2;
back   = + domain_size / 2;

left_inner   = 0.5*(inner_size*Cos(Pi/2 + alpha) - inner_size*Cos(alpha));
right_inner  = -left_inner;
back_inner    = 0.5*(inner_size*Sin(Pi/2 + alpha) + inner_size*Sin(alpha));
front_inner = -back_inner;

// outer domain boundary
Point(1) = {left,  back,    0, h};  // left top
Point(2) = {left,  front, 0, h};  // left front
Point(3) = {right, front, 0, h};  // right front
Point(4) = {right, back,    0, h};  // right top

// centers on boundaries
Point(5) = {left, 0, 0, h};  /* left */
Point(6) = {0, front, 0, h};  /* front */
Point(7) = {right, 0, 0, h};  /* right */
Point(8) = {0, back, 0, h};  /* back */

// inner domain boundary
Point(9) = {left_inner, 0.5*(inner_size*Sin(Pi/2 + alpha) - inner_size*Sin(alpha)), 0, hf}; // point on outer left
Point(10) = {left_inner + inner_size*Sin(alpha), front_inner, 0, hf}; // point on outer bottom
Point(11) = {right_inner, front_inner + inner_size*Sin(alpha), 0, hf}; // point on outer right
Point(12) = {left_inner + inner_size*Cos(alpha), back_inner, 0, hf}; // point on outer top

// inner domain internal points
Point(13) = {-inner_size/2*Cos(alpha), -inner_size/2*Sin(alpha), 0, hf};
Point(14) = {+inner_size/2*Cos(alpha), +inner_size/2*Sin(alpha), 0, hf};
Point(17) = {left_inner + (inner_size/2-lfrac/2)*Cos(alpha),
             0.5*(inner_size*Sin(Pi/2 + alpha) - inner_size*Sin(alpha)) + (inner_size/2-lfrac/2)*Sin(alpha), 0, hf};
Point(18) = {left_inner + (inner_size/2+lfrac/2)*Cos(alpha),
             0.5*(inner_size*Sin(Pi/2 + alpha) - inner_size*Sin(alpha)) +  (inner_size/2+lfrac/2)*Sin(alpha), 0, hf};
Point(19) = {left_inner + inner_size*Sin(alpha) + (inner_size/2-lfrac/2)*Cos(alpha),
             front_inner + (inner_size/2-lfrac/2)*Sin(alpha),
             0, hf}; // point on outer bottom
Point(20) = {left_inner + inner_size*Sin(alpha) + (inner_size/2+lfrac/2)*Cos(alpha),
             front_inner + (inner_size/2+lfrac/2)*Sin(alpha),
             0, hf}; // point on outer bottom

// frac
Point(15) = {-lfrac/2*Cos(alpha), -lfrac/2*Sin(alpha), 0, hf};
Point(16) = {+lfrac/2*Cos(alpha), +lfrac/2*Sin(alpha), 0, hf};

/* LINES */
/* outer domain */
Line(1) = {1, 5, 2}; // outer Left
Line(2) = {2, 6, 3}; // outer front
Line(3) = {3, 7, 4}; // outer right
Line(4) = {4, 8, 1}; // outer back

// inner domain
// upper half 1
Line(5) = {9, 13};
Line(6) = {13, 15};
Line(7) = {15, 16};  /* frac */
Line(8) = {16, 14};
Line(9) = {14, 12};
Line(10) = {12, 18};
Line(11) = {18, 17};
Line(12) = {17, 9};

Line(13) = {15, 17};
Line(14) = {18, 16};


// lower half
Line(15) = {13, 10};
Line(16) = {10, 19};
Line(17) = {19, 20};
Line(18) = {20, 11};
Line(19) = {11, 14};

Line(20) = {19, 15};
Line(21) = {16, 20};

// upper inner
Transfinite Line{5}  = (inner_size/2) / hf + 1;
Transfinite Line{6}  = (inner_size/2 - lfrac/2) / hf + 1;
Transfinite Line{7}  = nf + 1;
Transfinite Line{8}  = (inner_size/2 - lfrac/2) / hf + 1;
Transfinite Line{9}  = (inner_size/2) / hf + 1;
Transfinite Line{10}  = (inner_size/2 - lfrac/2) / hf + 1;
Transfinite Line{11}  = nf + 1;
Transfinite Line{12}  = (inner_size/2 - lfrac/2) / hf + 1;
Transfinite Line{13}  = (inner_size/2) / hf + 1;
Transfinite Line{14}  = (inner_size/2) / hf + 1;
/* Printf("%g ", (inner_size/2) / hf) ; */
/* Printf("%g ", (inner_size/2) / hf) ; */

// lower inner
Transfinite Line{15}  = (inner_size/2) / hf + 1;
Transfinite Line{16}  = (inner_size/2 - lfrac/2) / hf + 1;
Transfinite Line{17}  = nf + 1;
/* Transfinite Line{18}  = (inner_size/2 - lfrac/2) / hf + 1; */
Transfinite Line{19}  = (inner_size/2) / hf + 1;
Transfinite Line{20}  = (inner_size/2) / hf + 1;
Transfinite Line{21}  = (inner_size/2) / hf + 1;

/* // LOOPS */
Line Loop(1) = {1, 2, 3, 4};    // outer
// upper
Line Loop(2) = {5, 6, 13, 12};    // upper 1
Line Loop(3) = {-13, 7, -14, 11};    // upper 2
Line Loop(4) = {14, 8, 9, 10};    // upper 3
// lower
Line Loop(5) = {15, 16, 20, -6};    // upper 1
Line Loop(6) = {-20, 17, -21, -7};    // upper 1
Line Loop(7) = {21, 18, 19, -8};    // upper 1

Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7};      // outer domain
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};

// upper inner
Transfinite Surface(2) = {9, 13, 15, 17};
Transfinite Surface(3) = {17, 15, 16, 18};
Transfinite Surface(4) = {16, 14, 12, 18};

// lower inner
Transfinite Surface(5) = {13, 10, 19, 15};
Transfinite Surface(6) = {15, 19, 20, 16};
Transfinite Surface(7) = {16, 20, 11, 14};

// outer
Recombine Surface(1);
// upper inner
Recombine Surface(2);
Recombine Surface(3);
Recombine Surface(4);

// lower inner
Recombine Surface(5);
Recombine Surface(6);
Recombine Surface(7);

out1[] = Extrude{0,0,thickness} { Surface{1}; Layers{nz}; Recombine;};
out2[] = Extrude{0,0,thickness} { Surface{2}; Layers{nz}; Recombine;};
out3[] = Extrude{0,0,thickness} { Surface{3}; Layers{nz}; Recombine;};
out4[] = Extrude{0,0,thickness} { Surface{4}; Layers{nz}; Recombine;};
out5[] = Extrude{0,0,thickness} { Surface{5}; Layers{nz}; Recombine;};
out6[] = Extrude{0,0,thickness} { Surface{6}; Layers{nz}; Recombine;};
out7[] = Extrude{0,0,thickness} { Surface{7}; Layers{nz}; Recombine;};

// Labels
Left = 1111111;
Right = 1111112;
Front = 2222221;
Back = 2222222;
Bottom = 3333331;
Top = 3333332;

Frac1 = 1;

// Physical Surface(Bottom) = {1, 2, 3, 4, 5, 6, 7};
// Physical Surface(Top) = {163, 251, 185, 207, 229, 295, 273, 251};
// Physical Surface(Left) = {54};
// Physical Surface(Right) = {62};
// Physical Surface(Front) = {58};
// Physical Surface(Back) = {66};

Physical Surface(Frac1) = {94};

ELASTIC_DRIVER = 1;
SDA_DRIVER = 2;

Physical Volume(SDA_DRIVER) = {1, 2, 3, 4, 5, 6, 7};

Mesh.Algorithm = 8;

Geometry.Tolerance = size_x / nx / 1e3;
Geometry.AutoCoherence = 2;
// Mesh.MshFileVersion = 2.2;
Mesh.MshFileVersion = 4.1;

Mesh 3;  // Generalte 3D mesh
// Mesh.MshFileVersion = 4.1;
Mesh.MshFileVersion = 2.2;
Coherence Mesh;  // Remove duplicate entities
Save "geom.msh";  // Save mesh in MSH format
