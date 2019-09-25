#pragma once

#include "ConnectionMap.hpp"

#ifdef WITH_METIS
#include "metis.h"
#endif // WITH_METIS

#include <vector>
#include <exception>
#include <unordered_set>

namespace multiscale
{

using std::vector;
using std::size_t;

template <typename T>
using ConnectionMap = hash_algorithms::ConnectionMap<T>;


/* Wrapper class for METIS library */
template <typename ConnType>
class MetisInterface
{
 public:
  /* connections contains cell connections (or a graph)
   * n_blocks is how many partitions (coarse blocks) you wanna get
   * n_cells is how many unique untries are in the connection map;
   * if set to 0, then it's computed automatically */
  static vector<size_t> build_partitioning(const ConnectionMap<ConnType> & connections,
                                           const size_t n_blocks,
                                           size_t n_cells = 0)
  {
#ifdef WITH_METIS
    if (n_blocks < 2) throw std::invalid_argument("number of blocks too small");
    if (n_cells == 0) n_cells = count_elements(connections);

    // generate input: xadj (similar to row_ptr in CSR format)
    vector<idx_t> xadj(n_cells + 1, 0);
    for (auto it = connections.begin(); it != connections.end(); ++it)
    {
        const auto elements = it.elements();
        const int ia = elements.first;
        const int ib = elements.second;
        ++xadj[ia+1];
        ++xadj[ib+1];
    }

    // convert to the same ascending format as row_ptr
    for(std::size_t ib = 0; ib < n_cells; ++ib)
      xadj[ib+1] += xadj[ib];

    // generate input: adj (similar to col_ind in CSR format)
    vector<idx_t> adj(xadj[n_cells]);
    // graph connection weights
    vector<idx_t> adj_weight(xadj[n_cells]);
    for (auto it = connections.begin(); it != connections.end(); ++it)
    {
        const auto elements = it.elements();
        const int ia = elements.first;
        const int ib = elements.second;
        adj_weight[xadj[ia]] = small_wgt; //small weight between blocks
        adj[xadj[ia]++] = ib;
        adj_weight[xadj[ib]] = small_wgt; //small weight between blocks
        // temporarily change the value of xadj for generating adj
        adj[xadj[ib]++] = ia;
    }

    //  restore xadj
    for(std::size_t ib = n_cells; ib != 0; --ib)
      xadj[ib] = xadj[ib-1];
    xadj[0] = 0;

    // partition the entire domain: see the user manual of METIS for details
    idx_t ncon = 1, objval = 0;
    std::vector<idx_t> vwgt(n_cells, 1);
    std::vector<idx_t> size(n_cells, 1);
    idx_t icount = static_cast<idx_t>(n_cells);
    idx_t n_domains = static_cast<idx_t>(n_blocks);
    std::vector<real_t> ubvec(ncon, 1.1);

    // output: the corresponding thread of each grid block (default: 0th thread)
    vector<idx_t> coarse_cell_id(n_cells, 0);

    // call METIS graph partitioner
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NCUTS] = 100;
    METIS_PartGraphKway(&icount, // number of vertices in the graph
                        &ncon,   // n balalncing constraints (>= 1)
                        &xadj[0], &adj[0], // adjacency structure
                        &vwgt[0], &size[0],
                        NULL/*&adjwgt[0]*/,
                        &n_domains,
                        /* tpwgts = */NULL,  // weight for each partition and constraint
                        /* ubvec = */ &ubvec[0],
                        // /* ubvec = */NULL, // load imbalance tolerance
                        options, &objval,
                        &coarse_cell_id[0]);

    vector<size_t> partitioning(n_cells);
    for (std::size_t i=0; i<n_cells; ++i)
      partitioning[i] = static_cast<std::size_t>(coarse_cell_id[i]);

    return partitioning;

#else
  throw std::invalid_argument("METIS is not available, cannot perform partitioning");
#endif // WITH_METIS
  }  // end interface

 private:
  MetisInterface();

  // computes number of unique elements in connection map
  static size_t count_elements(const ConnectionMap<ConnType> & connections)
  {
    std::unordered_set<size_t> uniques;
    for (auto it = connections.begin(); it != connections.end(); ++it)
    {
      const auto elements = it.elements();
      uniques.insert(elements.first);
      uniques.insert(elements.second);
    }
    return uniques.size();
  }

  static const int small_wgt = 1;
  static const int large_wgt = 1000000;
};

}
