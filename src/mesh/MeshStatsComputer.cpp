#include "MeshStatsComputer.hpp"

namespace mesh {

MeshStatsComputer::MeshStatsComputer(const mesh::Mesh & grid)
    : _grid(grid)
{
  compute_edge_charachteristics_();
  compute_totals_();
}

void MeshStatsComputer::compute_totals_()
{
  _n_active_cells = _grid.n_active_cells();
  _n_cells = _grid.n_cells_total();
  _n_inactive_cells = _n_cells - _n_active_cells;
  _n_faces = _grid.n_faces();
  _n_inactive_faces = _n_faces - _n_active_faces;
}

void MeshStatsComputer::compute_edge_charachteristics_()
{
  // TODO: compute each edge only once
  _n_active_faces = 0;
  _total_edge_length = 0;
  _average_edge_length = 0;
  _maximum_edge_length = 0;
  _minimum_edge_length = std::numeric_limits<double>::max();
  size_t n_edges = 0;
  for (auto face = _grid.begin_active_faces(); face != _grid.end_active_faces(); ++face, ++_n_active_faces)
  {
    const auto & vertices = face->vertices();
    for (size_t i=0; i<vertices.size(); ++i)
    {
      const size_t i1 = vertices[i];
      const size_t i2 = (i == vertices.size() - 1) ? vertices[0] : vertices[i+1];
      const double l = _grid.vertex(i1).distance(_grid.vertex(i2));
      _total_edge_length += l;
      _maximum_edge_length = std::max(_maximum_edge_length, l);
      _minimum_edge_length = std::min(_minimum_edge_length, l);
      _average_edge_length += l;
      n_edges++;
    }
    _average_edge_length /= n_edges;
  }
}


}  // end namespace mesh
