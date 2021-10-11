#include "DiscretizationINSIM.hpp"

namespace discretization {

DiscretizationINSIM::DiscretizationINSIM(DoFNumbering const & dof_numbering,
                                         gprs_data::SimData & data,
                                         std::vector<ControlVolumeData> & cv_data,
                                         std::vector<ConnectionData> & connection_data)
    : DiscretizationBase(dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationINSIM::build()
{
  for (size_t v = 0; v < m_grid.n_vertices(); ++v)
    if ( m_dofs.has_vertex( v ) )
      build_vertex_data_( v );

  auto const vertex_adjacency = build_vertex_adjacency_();
  build_connectivity_( vertex_adjacency );
}

void DiscretizationINSIM::build_connectivity_(algorithms::EdgeWeightedGraph const & vertex_adjacency)
{
  for (size_t v = 0; v < m_grid.n_vertices(); ++v)
    for (auto const * e : vertex_adjacency.adj( v )) {
      size_t const w = e->other(v);
      if (v < w && m_dofs.has_vertex(v))
        std::cout << m_dofs.vertex_dof(v) << " " << m_dofs.vertex_dof(w) << std::endl;
    }
}

void DiscretizationINSIM::build_vertex_data_(size_t vertex)
{
  auto const attached_cells = m_grid.vertex_cells( vertex );
  size_t const n = attached_cells.size();

  std::vector<ControlVolumeData> attached_cvs( n );
  for (size_t i = 0; i < n; ++i)
    DiscretizationBase::build_cell_data_( *attached_cells[i], attached_cvs[i] );

  size_t const dof = m_dofs.vertex_dof( vertex );

  auto & cv = m_cv_data[dof];
  cv.type = ControlVolumeType::cell;
  cv.master = vertex;
  cv.center = m_grid.vertex( vertex );

  // take average of physical properties
  cv.custom.resize( attached_cvs.front().custom.size(), 0.f );
  for (size_t i = 0; i < n; ++i) {
    cv.volume += attached_cvs[i].volume;
    cv.porosity += attached_cvs[i].porosity / n;
    cv.permeability += (attached_cvs[i].permeability / (double)n);
    for (size_t j = 0; j < cv.custom.size(); ++j)
      cv.custom[j] += attached_cvs[i].custom[j] / n;
  }
}

algorithms::EdgeWeightedGraph DiscretizationINSIM::build_vertex_adjacency_() const
{
  algorithms::EdgeWeightedGraph g( m_grid.n_vertices() );
  for (auto cell = m_grid.begin_active_cells(); cell != m_grid.end_active_cells(); ++cell)
  {
    for (auto const & edge : cell->edges())
      if ( !g.has_edge( edge.first, edge.second ) )
        g.add( algorithms::UndirectedEdge(edge.first, edge.second, 0) );
  }
  return g;
}

}  // end namespace discretization
