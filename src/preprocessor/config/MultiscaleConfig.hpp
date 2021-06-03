#pragma once
#include <cstddef>  // size_t

enum class MSPartitioning : int
{
  no_partitioning  = 0,
  msrsb            = 1,  // igor's inspired by olav's paper, doesn't work for mech
  mrst_flow        = 2,   // jacques' inspired by mrst and cgal
  graph            = 3,   // igor's based on purely flow graph
  method_mechanics = 4  // igor's mechanics method
};

struct MultiscaleConfig
{
  MSPartitioning type = MSPartitioning::no_partitioning;
  size_t n_blocks = 2;
};
