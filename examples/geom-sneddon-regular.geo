// cubic domain with dimensions size_x, size_y, size_z
// number of cells nx, ny, nz
// make sure that ny is odd that it contains the central points
// of front and back surfaces in order to constrain
// y displacement

// pi = 3.1415926;
// nx = 50;
// ny = 51;
// nz = 1;

// outer domain size
size_x = 15;
size_y = 30;
lfrac = 2;

nf = 20;
h = lfrac / nf;
size_z = h;

nx = size_x / h;
ny = size_y / h;
nz = 1;
// make nfy odd
If(Floor(ny/2) == ny/2)
  ny = ny + 1;
EndIf

// Printf("%f", nfx/nfy);


// outer bomain points
left  = 0;
right = size_x;
back  = size_y/2;
front = -size_y/2;
bottom = 0;
top = size_z;

// domain boundary
Point(1) = {left,  back,  bottom, h};  // left top
Point(2) = {left,  front, bottom, h};  // left bottom
Point(3) = {right, front, bottom, h};  // right bottom
Point(4) = {right, back,  bottom, h};  // right top

// // LINES
Line(1) = {1, 2};    // left
Line(2) = {2, 3};    // front
Line(3) = {3, 4};    // right
Line(4) = {4, 1};    // back

// make frac cells structured
Transfinite Line{1} = ny + 1;
Transfinite Line{2} = nx + 1;
Transfinite Line{3} = ny + 1;
Transfinite Line{4} = nx + 1;

// LOOPS
Line Loop(1) = {1, 2, 3, 4};    // outer boundary

Plane Surface(1) = {1};      // outer surface

Transfinite Surface(1) = {1, 2, 3, 4}; // point labels

// remesh into quadrelaterals
Recombine Surface{1};

out1[] = Extrude{0,0,size_z} { Surface{1}; Layers{nz}; Recombine; };

// print
NumPoints = #out1[];
Printf("The Array Contents are") ;
For index In {0:NumPoints-1}
  Printf("%g ",out1[index]) ;
EndFor

// Labels
Left = 1111111;
Right = 1111112;
Front = 2222221;
Back = 2222222;
Bottom = 3333331;
Top = 3333332;

Physical Surface(Bottom) = {1};
Physical Surface(Top) = {out1[0]};
Physical Surface(Left) = {out1[2]};
Physical Surface(Right) = {out1[4]};
Physical Surface(Front) = {out1[3]};
Physical Surface(Back) = {out1[5]};

ELASTIC_DRIVER = 9999991;
SDA_DRIVER = 9999992;

Physical Volume(SDA_DRIVER) = {1};
