#pragma once

#include "angem/Point.hpp"
#include "angem/Tensor2.hpp"
#include <unordered_map>

namespace discretization {

enum ControlVolumeType
{
  cell,
  face
};

struct ControlVolumeData
{
  ControlVolumeType type;
  std::size_t       master; // grid cell/face index
  angem::Point<3,double> center;
  double volume = 0;
  double porosity = 0;
  double aperture = 0;  // faces only : fracture apreture
  angem::Tensor2<3,double> permeability;
  std::vector<double> custom;
};

}  // end namespace discretization
