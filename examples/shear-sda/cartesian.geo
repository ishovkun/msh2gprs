// cubic domain with dimensions size_x, size_y, size_z
// number of cells nx, ny, nz
// make sure that ny is odd that it contains the central points
// of front and back surfaces in order to constrain
// y displacement


// outer domain size
size_x = 60;
size_y = 60;

nx = 50;
ny = nx;
// nf = 4;
hx = size_x / nx;
hy = size_y / ny;
size_z = 1;

nz = 1;

// outer bomain points
left   = -size_x / 2;
right  = size_x  / 2;
back   = size_y  / 2;
front  = -size_y / 2;
bottom = -size_z / 2;
top    = size_z  / 2;

left_aq  = left  + hx;
right_aq = right - hx;

// domain boundary
Point(1) = {left,  back,  bottom, 1};  // left top
Point(2) = {left,  front, bottom, 1};  // left bottom
Point(3) = {right, front, bottom, 1};  // right bottom
Point(4) = {right, back,  bottom, 1};  // right top

Point(5) = {left_aq,  back,  bottom, 1};  // left top
Point(6) = {left_aq,  front, bottom, 1};  // left bottom
Point(7) = {right_aq, front, bottom, 1};  // right bottom
Point(8) = {right_aq, back,  bottom, 1};  // right top

// LINES
// left aquifer
Line(1) = {1, 2};    // left
Line(2) = {2, 6};    // front
Line(3) = {6, 5};    // right
Line(4) = {5, 1};    // back

// right aquifer
Line(5) = {8, 7};    // left
Line(6) = {7, 3};    // front
Line(7) = {3, 4};    // right
Line(8) = {4, 8};    // back

// main reservoir
Line(9) = {6, 7};    // front
Line(10) = {8, 5};    // back


// make frac cells structured
// left aq
Transfinite Line{1} = ny + 1;
Transfinite Line{2} = 1;
Transfinite Line{3} = ny + 1;
Transfinite Line{4} = 1;

// right aq
Transfinite Line{5} = ny + 1;
Transfinite Line{6} = 1;
Transfinite Line{7} = ny + 1;
Transfinite Line{8} = 1;

// main res
Transfinite Line{9} = nx + 1 - 2;
Transfinite Line{10} = nx + 1 - 2;

// LOOPS
Line Loop(1) = {1,  2, 3,  4};    // left aq
Line Loop(2) = {5,  6, 7,  8};    // right aq
Line Loop(3) = {-3, 9, -5, 10}; // main res

Plane Surface(1) = {1};      // left aq
Plane Surface(2) = {2};      // right aq
Plane Surface(3) = {3};      // main res

Transfinite Surface(1) = {1, 2, 6, 5}; // left aq
Transfinite Surface(2) = {8, 7, 3, 4}; // right aq
Transfinite Surface(3) = {5, 6, 7, 8}; // right aq

/* // remesh into quadrelaterals */
Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};

out1[] = Extrude{0,0,size_z} { Surface{1}; Layers{nz}; Recombine; };
out2[] = Extrude{0,0,size_z} { Surface{2}; Layers{nz}; Recombine; };
out3[] = Extrude{0,0,size_z} { Surface{3}; Layers{nz}; Recombine; };

// print
/* NumPoints = #out1[]; */
/* Printf("The Array Contents are") ; */
/* For index In {0:NumPoints-1} */
/*   Printf("%g ",out1[index]) ; */
/* EndFor */

// Labels
Left = 1111111;
Right = 1111112;
Front = 2222221;
Back = 2222222;
Bottom = 3333331;
Top = 3333332;

Physical Surface(Bottom) = {1, 2, 3};
Physical Surface(Top) = {32, 76, 54};
Physical Surface(Left) = {19};
Physical Surface(Right) = {49};
Physical Surface(Front) = {23, 67, 45};
Physical Surface(Back) = {31, 75, 53};

ELASTIC_DRIVER = 9999991;
SDA_DRIVER = 9999992;
AQUIFER = 9999993;
AQUIFER2 = 9999994;

Physical Volume(SDA_DRIVER) = {1, 2, 3};
