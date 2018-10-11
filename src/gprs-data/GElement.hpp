#pragma once

#include <Point.hpp>

struct Gelement
{
  int nMarker;
  int nVertices;
  int nNeighbors;

  int vtkIndex, formIndex;
  int fluidElement;

  std::vector<std::size_t> vVertices;
  std::vector<std::size_t> vVerticesSorted;
  std::vector<int> vVerticesNewnum;

  std::vector<std::size_t> vNeighbors;

  double thickness, center_distance;
  angem::Point<3,double> center;
  angem::Point<3,double> normal;
  double aperture, conductivity;
};
