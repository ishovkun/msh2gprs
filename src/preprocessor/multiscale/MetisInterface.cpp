#include "MetisInterface.hpp"

namespace multiscale {

std::vector<size_t>
MetisInterface::partition(algorithms::EdgeWeightedGraph const &g, size_t n_blocks)
{
  if (n_blocks < 2) throw std::invalid_argument("number of blocks too small");

  // generate input: xadj (similar to row_ptr in CSR format)
  std::vector<idx_t> xadj(g.n_vertices() + 1, 0);
  for (auto const &edge : g.edges()) {
    size_t const v = edge.either();
    size_t const w = edge.other(v);
    ++xadj[v + 1];
    ++xadj[w + 1];
  }
  // convert to the same ascending format as row_ptr
  for (std::size_t ib = 0; ib < g.n_vertices(); ++ib)
    xadj[ib + 1] += xadj[ib];

  // generate input: adj (similar to col_ind in CSR format)
  vector<idx_t> adj(xadj[g.n_vertices()]);
  // graph connection weights
  vector<idx_t> adj_weight(xadj[g.n_vertices()]);
  for (auto const &edge : g.edges()) {
    size_t const v = edge.either();
    size_t const w = edge.other(v);
    adj_weight[xadj[v]] = edge.weight();
    adj[xadj[v]++] = w;
    adj_weight[xadj[w]] = edge.weight();
    adj[xadj[w]++] = v;
  }

  //  restore xadj
  for (size_t v = g.n_vertices(); v > 0; --v)
    xadj[v] = xadj[v - 1];
  xadj[0] = 0;

  // output: the corresponding thread of each grid block (default: 0th thread)
  vector<idx_t> coarse_cell_id(g.n_vertices(), 0);
#ifdef WITH_METIS
  // partition the entire domain: see the user manual of METIS for details
  idx_t ncon = 1, objval = 0;
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NCUTS] = 8;
  std::vector<idx_t> vwgt(g.n_vertices(), 1);
  std::vector<idx_t> size(g.n_vertices(), 1);
  idx_t icount = static_cast<idx_t>(g.n_vertices());
  idx_t n_domains = static_cast<idx_t>(n_blocks);

  // call METIS graph partitioner
  METIS_PartGraphKway(&icount, &ncon, &xadj[0], &adj[0], &vwgt[0], &size[0],
                      NULL /*&adjwgt[0]*/, &n_domains, NULL, NULL, options,
                      &objval, &coarse_cell_id[0]);
#else
  throw std::invalid_argument(
      "METIS is not available, cannot perform partitioning");
#endif // WITH_METIS

  vector<size_t> part(g.n_vertices());
  for (size_t i = 0; i < g.n_vertices(); ++i)
    part[i] = static_cast<std::size_t>(coarse_cell_id[i]);

  return part;
}

}  // end namespace multiscale
