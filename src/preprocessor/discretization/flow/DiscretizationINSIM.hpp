#pragma once
#include "DiscretizationBase.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"

namespace discretization {

class DiscretizationINSIM : public DiscretizationBase {
 public:
  // constructor
  DiscretizationINSIM(DoFNumbering const & dof_numbering,
                     gprs_data::SimData & data,
                     std::vector<ControlVolumeData> & cv_data,
                     std::vector<ConnectionData> & connection_data);

  // main method. build the discretization
  void build() override;

  // destructor
  virtual ~DiscretizationINSIM() = default;

 private:
  void build_vertex_data_(size_t vertex);
  algorithms::EdgeWeightedGraph build_vertex_adjacency_() const;
  // build vector of connections
  void build_dof_adjecency_(algorithms::EdgeWeightedGraph const & vertex_adjacency);
  std::vector<size_t> bfs_(size_t u, algorithms::EdgeWeightedGraph const & vertex_adjacency);
};

}  // end namespace discretization
