// pi = 3.1415926;
alpha = 30 * Pi / 180;

domain_size = 200; // outer
l = 10;

inner_size = 5*l;
thickness = -10;

nf = 4;
h = l / nf;  // outer boundary element size
hf = l / nf;   // fracture element size - width

left  = -domain_size / 2;
right = +domain_size / 2;
front = -domain_size / 2;
back   = +domain_size / 2;

left_inner   = 0.5*(inner_size*Cos(Pi/2 + alpha) - inner_size*Cos(alpha));
right_inner  = -left_inner;
back_inner    = 0.5*(inner_size*Sin(Pi/2 + alpha) + inner_size*Sin(alpha));
front_inner = -back_inner;

// outer domain boundary
Point(1) = {left,  back,    0, h};  // left top
Point(2) = {left,  front, 0, h};  // left front
Point(3) = {right, front, 0, h};  // right front
Point(4) = {right, back,    0, h};  // right top


// inner domain boundary
Point(5) = {left_inner, 0.5*(inner_size*Sin(Pi/2 + alpha) - inner_size*Sin(alpha)), 0, hf}; // point on outer left
Point(6) = {left_inner + inner_size*Sin(alpha), front_inner, 0, hf}; // point on outer bottom
Point(7) = {right_inner, front_inner + inner_size*Sin(alpha), 0, hf}; // point on outer right
Point(8) = {left_inner + inner_size*Cos(alpha), back_inner, 0, hf}; // point on outer top

// LINES
// outer domain
Line(1) = {1, 2}; // outer Left
Line(2) = {2, 3}; // outer front
Line(3) = {3, 4}; // outer right
Line(4) = {4, 1}; // outer back

// inner domain
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

n_perp = Round(inner_size / hf);
n_along = Round(inner_size / hf);
If (n_perp % 2 == 0)
  n_perp += 1;
EndIf
If (n_along % 2 != 0)
    n_along += 1;
EndIf
Printf("%d", n_along);
Printf("%d", n_perp);

// make frac cells structured
Transfinite Line{5}  = n_perp + 1;
Transfinite Line{6}  = n_along+ 1;
Transfinite Line{7}  = n_perp + 1;
Transfinite Line{8}  = n_along + 1;

// LOOPS
Line Loop(1) = {1, 2, 3, 4};    // outer
Line Loop(2) = {5, 6, 7, 8};    // inner

Plane Surface(1) = {1, 2};      // outer domain left top
Plane Surface(2) = {2};      // outer domain left top
Transfinite Surface(2) = {5, 6, 7, 8};

Recombine Surface(1);
Recombine Surface(2);

Color Green{ Surface{ 1 }; }
Color Purple{ Surface{ 2 }; }

out1[] = Extrude{0,0,thickness} { Surface{1}; Layers{1}; Recombine;};
out1[] = Extrude{0,0,thickness} { Surface{2}; Layers{1}; Recombine;};

// Labels
Left = 1111111;
Right = 1111112;
Front = 2222221;
Back = 2222222;
Bottom = 3333331;
Top = 3333332;

Physical Surface(Bottom) = {1, 2};
Physical Surface(Top) = {50, 72};
Physical Surface(Left) = {21};
Physical Surface(Right) = {29};
Physical Surface(Front) = {25};
Physical Surface(Back) = {33};

/* ELASTIC_DRIVER = 9999991; */
SDA_DRIVER = 9999992;

Physical Volume(SDA_DRIVER) = {1, 2};

// For even better quadrilateral meshes, you can try the experimental "Delaunay
// for quads" (DelQuad) meshing algorithm: DelQuad is a triangulation algorithm
// that enables to create right triangles almost everywhere. Uncomment the
// following line to try DelQuad:
Mesh.Algorithm = 8;
