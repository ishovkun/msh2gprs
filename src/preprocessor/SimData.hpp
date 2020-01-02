#pragma once

#include "mesh/Mesh.hpp"
#include "angem/Tensor2.hpp"

namespace gprs_data {

struct SimData
{
  mesh::Mesh grid;
  // cell properties
  std::vector<std::string> property_names;
  std::vector<std::vector<double>> cell_properties;
  // permeability keys
  std::array<double,9> permeability_keys;
  size_t porosity_key_index;
  std::vector<size_t> output_flow_properties;


  angem::Tensor2<3,double> get_permeability(const std::size_t cell) const
  {
    assert(cell < cell_properties[permeability_keys[0]].size());
    angem::Tensor2<3,double> K;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        K(i, j) = (permeability_keys[3*i+j] >= 0) ?
                  cell_properties[permeability_keys[3*i+j]][cell] : 0;
    return K;
  }

  double get_porosity(const std::size_t cell) const
  {
    assert(cell < property_names.size());
    assert(porosity_key_index > 0 && porosity_key_index < property_names.size());
    return cell_properties[porosity_key_index][cell];
  }
};


}  // end namespace gprs_data
