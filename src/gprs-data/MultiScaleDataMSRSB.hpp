#pragma once

#include "mesh/Mesh.hpp"
#include "LayerDataMSRSB.hpp"
#include <algorithm>  // std::max_element
#include <vector>

namespace multiscale
{

using std::size_t;
using std::vector;

class MultiScaleDataMSRSB
{
 public:
  /* Constructor.
   * takes n_blocks for only a single layer,
   * since multi-level multiscale is a long way
   * down the road. */
  MultiScaleDataMSRSB(mesh::Mesh  & grid,
                      const size_t  n_blocks);
  // get reference to the active layer
  LayerDataMSRSB & active_layer(){return layers[active_layer_index];}

 protected:
  void build_partitioning();

 private:
  const mesh::Mesh & grid;
  vector<LayerDataMSRSB> layers;
  size_t active_layer_index;
};


// void MultiScaleDataMSRSB::()
// {

// }


}
