#include "DiscretizationINSIM.hpp"
#include <queue>

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
  build_dof_adjecency_( vertex_adjacency );
}

void DiscretizationINSIM::build_dof_adjecency_(algorithms::EdgeWeightedGraph const & vertex_adjacency)
{
  algorithms::EdgeWeightedGraph g( m_dofs.n_dofs() );


  for (size_t u = 0; u < m_grid.n_vertices(); ++u)
    if ( m_dofs.has_vertex( u ) ) {
      std::vector<size_t> neighbors = bfs_(u, vertex_adjacency);
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

struct QueueItem {
  size_t vertex;
  size_t level;
};

std::vector<size_t> DiscretizationINSIM::bfs_(size_t source, algorithms::EdgeWeightedGraph const & vertex_adjacency)
{
  std::queue<QueueItem> q;
  std::vector<bool> visited( vertex_adjacency.n_vertices(), false );
  q.push({source, 0});
  visited[source] = true;
  std::vector<size_t> neighbors;

  size_t max_level = std::numeric_limits<size_t>::max();
  bool found = false;
  size_t cnt = 0;
  while ( !q.empty() && q.front() < max_level ) {
    // take item from the queue
    auto const item = q.front();
    size_t const u = item.vertex;
    q.pop();

    // check if we found a correct dof to connect to
    if ( m_dofs.has_vertex( u ) && u != source ) {
      found = true;
      neighbors.push_back( u );
      max_level = item.level;  // don't search beyond this depth
    }
    else
    {
      for (auto * edge : vertex_adjacency.adj(u)) {  // search neighbors
        size_t const v = edge->other(u);
        if ( !visited[v] ) {
          visited[v] = true;
          q.push({ v,  item.level + 1 });
        }
      }
    }
  }
}

}  // end namespace discretization
