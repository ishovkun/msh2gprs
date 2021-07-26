#pragma once
#include <array>
#include <cstddef>  // size_t

enum class MSPartitioning : int
{
  no_partitioning = 0,
  metis           = 1,
  geometric       = 2,
};

enum class MSSupportType
{
  msrsb     = 1, // igor's inspired by olav's paper, doesn't work for mech
  graph     = 3, // igor's based on purely flow graph
  mechanics = 4         // igor's mechanics method
};

struct MultiscaleConfig
{
  MSPartitioning part_type = MSPartitioning::no_partitioning;
  MSSupportType support_type = MSSupportType::graph;
  // if geometric, we need three values
  // if metis, we just take the first value
  std::array<size_t,3> n_blocks{2, 1, 1};
};
