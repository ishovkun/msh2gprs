#include "DiscretizationINSIM.hpp"
#include <queue>

#include "mesh/io/VTKWriter.hpp"

namespace discretization {

using namespace algorithms;

DiscretizationINSIM::DiscretizationINSIM(DoFNumbering const & dof_numbering,
                                         DoFNumbering const & vertex_to_well,
                                         gprs_data::SimData & data,
                                         std::vector<ControlVolumeData> & cv_data,
                                         std::vector<ConnectionData> & connection_data)
    : DiscretizationBase(dof_numbering, data, cv_data, connection_data)
    , _vertex_to_well(vertex_to_well)
{}

void DiscretizationINSIM::build()
{
  for (size_t v = 0; v < m_grid.n_vertices(); ++v) {
    if ( m_dofs.has_vertex( v ) ) {
      build_vertex_data_( v );
    }
  }

  auto const vertex_adjacency = build_vertex_adjacency_();
  auto const dof_adjacency = build_dof_adjecency_( vertex_adjacency );
  build_connections_(dof_adjacency);
}

void DiscretizationINSIM::build_connections_(algorithms::EdgeWeightedGraph const & dof_adjacency)
{
  for (auto const & edge : dof_adjacency.edges()) {
    // we really need tube volume, transmissibility, and tube length
    m_con_data.emplace_back();
    auto & con = m_con_data.back();
    con.type = ConnectionType::matrix_matrix;
    size_t const dof_u = edge.either();
    size_t const dof_v = edge.other(dof_u);
    con.elements = {dof_u, dof_v};
    auto const & cell1 = m_cv_data[dof_u];
    auto const & cell2 = m_cv_data[dof_v];
    // std::cout << "building con " << cell1.master << "-" << cell2.master << std::endl;
    con.center = 0.5 * (cell1.center + cell2.center);
    // project permeability
    auto const & K1 = cell1.permeability;
    auto const & K2 = cell2.permeability;
    // directional permeability
    auto const & cp = con.center;
    auto const & c1 = cell1.center;
    auto const & c2 = cell2.center;
    const double Kp1 = (K1 * (c1 - cp).normalize()).norm();
    const double Kp2 = (K2 * (c2 - cp).normalize()).norm();

    double const poro = 0.5 * ( cell1.porosity + cell2.porosity );
    double const d = cell1.center.distance( cell2.center );
    // double const pore_volume = approximate_connection_pore_volume_(cell1, cell2);
    // double const connection_volume = pore_volume / poro;
    double const connection_volume = d * (d/2 * d/2);
    double const pore_volume = connection_volume * poro;
    con.area = connection_volume / d;  // approximately :-)
    con.distances = {d/2, d/2};        // storage for tube length

    const double T1 = con.area * Kp1 / con.distances[0];
    const double T2 = con.area * Kp2 / con.distances[1];

    // face transmissibility
    const double T = T1*T2 / ( T1 + T2 );
    con.coefficients = {T, pore_volume};
  }
}

algorithms::EdgeWeightedGraph DiscretizationINSIM::build_dof_adjecency_(algorithms::EdgeWeightedGraph const & vertex_adjacency)
{
  EdgeWeightedGraph g( m_dofs.n_dofs() );

  /*
  ** Start a bfs from each dof vertex.
  ** Search until we found a connection at a particular search depth "level".
  ** Even if we found a well, keep searching until we look at all wells at this depth, and then stop.
  ** NOTE: there is room for optimization, we probably should not consider paths with vertices connected previously.
   */
  for (size_t u = 0; u < m_grid.n_vertices(); ++u)
    if ( m_dofs.has_vertex( u ) ) {
      std::vector<size_t> const neighbors = bfs_(u, vertex_adjacency);
      size_t const dof_u = m_dofs.vertex_dof( u );
      for (size_t const v : neighbors) {
        size_t const dof_v = m_dofs.vertex_dof( v );
        if ( !g.has_edge( dof_u, dof_v ) )
          g.add( UndirectedEdge(dof_u, dof_v, 0.f) );
      }
    }

  return g;
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

  assert( !attached_cvs.empty() );
  cv.custom.resize( attached_cvs.front().custom.size(), 0.f );
  for (size_t i = 0; i < n; ++i) {
    double const vi = attached_cvs[i].volume;
    cv.volume += vi;
    cv.porosity += attached_cvs[i].porosity * vi;
    cv.permeability += attached_cvs[i].permeability * vi;
    for (size_t j = 0; j < cv.custom.size(); ++j)
      cv.custom[j] += attached_cvs[i].custom[j] * vi;
  }

  // take average (divide by tot volume)
  cv.porosity /= cv.volume;
  cv.permeability /= cv.volume;
    for (size_t j = 0; j < cv.custom.size(); ++j)
      cv.custom[j] /= cv.volume;

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
  visited[source] = true;
  q.push({ source, 0 });
  std::vector<size_t> neighbors;
  size_t const source_well = _vertex_to_well.vertex_dof( source );

  size_t max_level = std::numeric_limits<size_t>::max();
  bool found = false;
  size_t cnt = 0;

  while ( !q.empty() && q.front().level <= max_level ) {
    // take item from the queue
    auto const item = q.front();
    size_t const u = item.vertex;
    q.pop();

    // check if vertex is a flow dof,
    // additionally check that we don't connect vertices that belong to the same well
    bool const same_well = _vertex_to_well.has_vertex(u) && _vertex_to_well.vertex_dof(u) == source_well;
    bool const same_dof = source == u;
    if ( m_dofs.has_vertex( u ) && !same_well && !same_dof ) {
      found = true;
      neighbors.push_back( u );
      max_level = item.level;  // don't search beyond this depth
    }
    else
    {
      for (auto const * const edge : vertex_adjacency.adj(u)) {  // search neighbors
        size_t const v = edge->other(u);
        if ( !visited[v] ) {
          visited[v] = true;
          q.push({ v,  item.level + 1 });
        }
      }
    }
  }

  return neighbors;
}

double DiscretizationINSIM::approximate_connection_pore_volume_(ControlVolumeData const & u, ControlVolumeData const & v)
{
  angem::LineSegment<double> segment( u.center, v.center );
  double vol = 0.f;
  for (const size_t cell_index : m_data.grid_searcher->collision(segment)) {
    const auto & cell = m_data.grid.cell(cell_index);
    double poro = m_data.cell_properties[m_data.flow.porosity_idx][cell_index];
    vol += cell.volume() * poro;
  }

  assert( vol > 0 && "Bug in volume computation");
  return vol;
}


}  // end namespace discretization
