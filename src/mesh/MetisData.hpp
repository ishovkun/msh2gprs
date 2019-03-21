#pragma once

#include <vector>


namespace mesh
{

struct MetisData
{
  // The number of balancing constraints. It should be at least 1
  idx_t ncon;
  // The adjacency structure of the graph
  std::vector<idx_t> xadj, adj;
  // The weights of the vertices as described in Section 5.5.
  std::vector<idx_t> vwgt;
  // The size of the vertices for computing the total communication
  // volume as described in Section 5.7.
  std::vector<idx_t> size;
  // The weights of the edges as described in Section 5.5.
  std::vector<idx_t> adjwgt;
};

/*
  // Exmplanation of METIS parameters
  nvtxs
  The number of vertices in the graph.

  ncon
  The number of balancing constraints. It should be at least 1.

  xadj, adjncy
  The adjacency structure of the graph (Section 5.5 in)

  vwgt (NULL)
  The weights of the vertices as described in Section 5.5.

  vsize (NULL)
  The size of the vertices for computing the total communication
  volume as described in Section 5.7.

  adjwgt (NULL)
  The weights of the edges as described in Section 5.5.

  nparts
  The number of parts to partition the graph.

  tpwgts (NULL)
  This is an array of size nparts×ncon that specifies the
  desired weight for each partition and constraint.
  The target partition weight for the ith partition and jth constraint is
  specified at tpwgts[i*ncon+j] (the numbering for both partitions and
  constraints starts from 0). For each constraint, the sum of the
  tpwgts[] entries must be 1.0


  ubvec (NULL)
  This is an array of size ncon that specifies the allowed load
  imbalance tolerance  for each  constraint. For the  ith partition  and jth
  constraint the allowed weight is the ubvec[j]*tpwgts[i*ncon+j] fraction of
  the jth’s  constraint total  weight. The load  imbalances must  be greater
  than 1.0.  A NULL value can  be passed indicating that  the load imbalance
  tolerance for  each constraint should be  1.001 (for ncon=1) or  1.01 (for
  ncon¿1).

  options (NULL)
  This is the array of options as described in Section 5.4.
  The following options are valid for METIS PartGraphRecursive:
  METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE, METIS_OPTION_RTYPE,
  METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS, METIS_OPTION_NITER,
  METIS_OPTION_SEED, METIS_OPTION_UFACTOR, METIS_OPTION_NUMBERING,
  METIS_OPTION_DBGLVL

  objval
  Upon successful completion, this variable stores the edge-cut or the total communication volume of
  the partitioning solution. The value returned depends on the partitioning’s objective function.

  part
  This is a vector of size nvtxs that upon successful completion stores the partition vector of the graph.
  The numbering of this vector starts from either 0 or 1, depending on the value of
  options[METIS OPTION NUMBERING].
*/

}  // end namespace
