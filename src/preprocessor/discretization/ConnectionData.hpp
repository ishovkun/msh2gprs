#pragma once

#include "angem/Point.hpp"
#include <vector>

namespace discretization {

enum ConnectionType
{
  matrix_matrix = 1,
  matrix_fracture = 2,
  fracture_fracture = 3
};

struct ConnectionData
{
  ConnectionType type;
  std::vector<double> coefficients;  // transmissibilities
  std::vector<size_t> elements;      // dofs that form a connection
  std::vector<double> distances;     // distance from cv center to connection
  double area;                       // connection area
  angem::Point<3,double> normal;
  angem::Point<3,double> center;
  // for fractures only
  double edge_length;                              // frac-frac edge length
  // angem::Point<3,double> projection;               // projection of cell/face onto connecting face/edge
  // angem::Point<3,double> edge_direction;           // normalizer unit vector from one vertex to another
  std::vector<size_t> all_elements;      // star delta participants
};

}  // end namespace discretization
