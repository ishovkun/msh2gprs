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
  // map control volume to cell/face index
  std::size_t       master;
  angem::Point<3,double> center;
  double volume;
  double porosity;
  angem::Tensor2<3,double> permeability;
  std::vector<double> custom;
};

}  // end namespace discretization
