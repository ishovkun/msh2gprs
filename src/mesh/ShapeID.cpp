#include <ShapeID.hpp>

namespace mesh
{

ShapeIndexer::MapIntInt ShapeIndexer::vtk_indices = {
  {TRGLE3,  5},
  {QUAD4,   9},
  {TETRA4,  10},
  {PRISM8,  12},
  {PRISM6,  13},
  {TRGLE6,  22},
  {QUAD8,   23},
  {PRISM20, 25},
  {PRISM15, 26}
};


ShapeIndexer::MapIntInt ShapeIndexer::shape_ids = {
  {4, TETRA4},
  {6, PRISM6}
};

}
