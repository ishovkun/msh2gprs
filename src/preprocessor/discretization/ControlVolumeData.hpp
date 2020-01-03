#pragma once

#include "angem/Point.hpp"
#include "angem/Tensor2.hpp"

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
  double volume;
  double porosity;
  angem::Tensor2<3,double> permeability;
  std::vector<double> custom;
};

}  // end namespace discretization
