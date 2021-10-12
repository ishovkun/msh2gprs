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
  // build connection graph for all vertices of the grid
  algorithms::EdgeWeightedGraph build_vertex_adjacency_() const;
  // build connection graph for all dofs
  algorithms::EdgeWeightedGraph build_dof_adjecency_(algorithms::EdgeWeightedGraph const & vertex_adjacency);
  // get a list of connected dofs (grid vertex indices not really dofs) for vertex u
  std::vector<size_t> bfs_(size_t u, algorithms::EdgeWeightedGraph const & vertex_adjacency);
  // build std::vector<ConnectionData> for output
  void build_connections_(algorithms::EdgeWeightedGraph const & dof_adjacency);
  // compute which cells does u-v line segment intersect, and add up their pore volume
  double approximate_connection_pore_volume_(ControlVolumeData const & u, ControlVolumeData const & v);
};

}  // end namespace discretization
