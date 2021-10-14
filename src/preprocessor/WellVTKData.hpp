#pragma once
#include "angem/PointSet.hpp"

struct WellVTKGrid
{
  angem::PointSet<3,double> vertices;  // set of well coordinatees: used for vtk output.
  std::vector<std::pair<std::size_t,std::size_t>> indices;
};
