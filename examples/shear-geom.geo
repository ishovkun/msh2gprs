// pi = 3.1415926;
alpha = 30 * Pi / 180;
l = 3;
h = 4*l;  // outer boundary element size


inl = 1.5*2*l; // size of inner domain along the fracture
inw = 2*l; // size of inner domain perpendicular the fracture

// hf = l/10;   // fracture element size - width
nf = 1.5*16;
nf1 = Floor(nf * inw/inl);

hf = inl / nf;

If(Floor(nf1/2) == nf1/2)
  nf1 = nf1 + 1;
EndIf

left  = -30*l;
right = +30*l;
back  = +30*l;
front = -30*l;

// outer domain boundary
Point(1) = {left,  back,    0, h};  // left top
Point(2) = {left,  front, 0, h};  // left bottom
Point(3) = {right, front, 0, h};  // right bottom
Point(4) = {right, back,    0, h};  // right top

// central points on left and right boundaries (for y displacement constraints)
Point(5) = {left,  0,    0, h};  // left center
Point(6) = {right,  0,    0, h};  // right center

// // inner domain boundary
Point(7) = {+inl/2*Cos(alpha) + inw/2*Cos(Pi/2+alpha), +inl/2*Sin(alpha) + inw/2*Sin(Pi/2+alpha), 0, hf}; // top left
Point(8) = {-inl/2*Cos(alpha) + inw/2*Cos(Pi/2+alpha), -inl/2*Sin(alpha) + inw/2*Sin(Pi/2+alpha), 0, hf}; // bot left
Point(9) = {-inl/2*Cos(alpha) - inw/2*Cos(Pi/2+alpha), -inl/2*Sin(alpha) - inw/2*Sin(Pi/2+alpha), 0, hf}; // bot right
Point(10) = {+inl/2*Cos(alpha) - inw/2*Cos(Pi/2+alpha), +inl/2*Sin(alpha) - inw/2*Sin(Pi/2+alpha), 0, hf}; // top right

// LINES
// outer domain
Line(1) = {1, 5, 2}; // outer left
Line(2) = {2, 3};    // outer bottom
Line(3) = {3, 6, 4}; // outer right
Line(4) = {4, 1};    // outer bottom
// inner domain
Line(5) = {7, 8};     // outer left
Line(6) = {8, 9};     // outer bottom
Line(7) = {9, 10};    // outer right
Line(8) = {10, 7};    // outer top


// make frac cells structured
Transfinite Line{5} = nf + 1;
Transfinite Line{6} = nf1 + 1;
Transfinite Line{7} = nf + 1;
Transfinite Line{8} = nf1 + 1;

// LOOPS
Line Loop(1) = {1, 2, 3, 4};    // outer boundary
Line Loop(2) = {5, 6, 7, 8};    // inner boundary

Plane Surface(1) = {1, 2};      // outer surface
Plane Surface(2) = {2};      // inner surface

Transfinite Surface(2) = {7, 8, 9, 10}; // point labels

// remesh into quadrelaterals
Recombine Surface{2};
Recombine Surface{1};

// // // Color Green{ Surface{ 1 }; }
// // // Color Purple{ Surface{ 2 }; }

out1[] = Extrude{0,0,hf} { Surface{1}; Layers{1}; Recombine; };
out2[] = Extrude{0,0,hf} { Surface{2}; Layers{1}; Recombine; };

// Labels
Left = 1111111;
Right = 1111112;
Front = 2222221;
Back = 2222222;
Bottom = 3333331;
Top = 3333332;

Physical Surface(Bottom) = {1, 2};
Physical Surface(Top) = {out1[0], out2[0]};
Physical Surface(Left) = {out1[2]};
Physical Surface(Right) = {out1[4]};
Physical Surface(Front) = {out1[3]};
Physical Surface(Back) = {out1[5]};

ELASTIC_DRIVER = 9999991;
SDA_DRIVER = 9999992;

Physical Volume(SDA_DRIVER) = {1, 2};
