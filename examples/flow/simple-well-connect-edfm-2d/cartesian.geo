// cubic domain with dimensions size_x, size_y, size_z
// number of cells nx, ny, nz
// make sure that ny is odd that it contains the central points
// of front and back surfaces in order to constrain
// y displacement


// outer domain size
size_x = 100;
size_y = 100;
size_z = 10;
depth = 0;

nx = 10;
ny = 10;
nz = 2;

hx = size_x / nx;
hy = size_y / ny;

// outer bomain points
left   = -size_x / 2;
right  = size_x  / 2;
back   = size_y  / 2;
front  = -size_y / 2;
bottom = depth - size_z;
top    = depth;

// domain boundary
Point(1) = {left,  back,  bottom, 1};  // left top
Point(2) = {left,  front, bottom, 1};  // left bottom
Point(3) = {right, front, bottom, 1};  // right bottom
Point(4) = {right, back,  bottom, 1};  // right top

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

Flow_Domain = 1;

Physical Volume(Flow_Domain) = {1};

Mesh 3;  // Generalte 3D mesh
Mesh.MshFileVersion = 2.2;
// Mesh.Algorithm = 8;
// Mesh.Algorithm3D = 1;
Coherence Mesh;  // Remove duplicate entities
Save "grid.msh";  // Save mesh in MSH format
// // Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;
