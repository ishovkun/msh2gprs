#pragma once

#include <unordered_map>

namespace mesh
{

enum ShapeID
{
  TRGLE3  = 0,
  TRGLE6  = 1,
  QUAD4   = 2,
  TETRA4  = 3,
  TETRA10 = 4,
  PRISM6  = 5,
  PRISM8  = 6,
  QUAD8   = 7,
  PRISM15 = 8,
  PRISM20 = 9
};


class ShapeIndexer
{
 public:
  static int get_vtk_index(const int shape_index) {return vtk_indices[shape_index];}
  static int get_standard_shape (const int n_vertices) {return shape_ids[n_vertices];}

  typedef std::unordered_map<int,int> MapIntInt;
  // map shape_id -> vtk index
  static MapIntInt vtk_indices;
  // map n_vertices -> shape_id
  static MapIntInt shape_ids;
};



}
