#include "FractureTip.hpp"

namespace mesh {

FractureTip::FractureTip(const SurfaceMesh<double> & face_mesh,
                         const Mesh & grid,
                         std::unordered_map<size_t,size_t> & surface_vertex_to_global)
    : _face_mesh(face_mesh), _grid(grid), _surface_vertex_to_global(surface_vertex_to_global)
{
  _global_boundary = find_global_boundary_();
  _tip_verts.resize( _grid.n_vertices() );
  for (auto edge = face_mesh.begin_edges(); edge != face_mesh.end_edges(); ++edge)
    if (edge.neighbors().size() == 1)
    {
      const auto edge_vertices = edge.vertex_indices();
      const size_t v1 = _surface_vertex_to_global[edge_vertices.first];
      const size_t v2 = _surface_vertex_to_global[edge_vertices.second];
      _tip_verts[v1].push_back(v2);
      _tip_verts[v2].push_back(v1);
    }
}

std::unordered_set<size_t> FractureTip::find_global_boundary_()
{
  std::unordered_set<size_t> result;
  for (auto face = _grid.begin_active_faces(); face != _grid.end_active_faces(); ++face)
    if ( face->neighbors().size() == 1 )
      for (auto v : face->vertices())
        result.insert( v );
  return result;
}

bool FractureTip::contains(const size_t vertex)
{
  if (vertex >= _tip_verts.size())
    throw std::invalid_argument("vertex " + std::to_string(vertex) + " does not exist");

  if (_tip_verts[vertex].empty())
    return false;
  else
  {
    if (_global_boundary.find(vertex) == _global_boundary.end())
      return true;
    else
    {
      // for (const size_t jvertex : _tip_verts[vertex])
      //   if ( _global_boundary.find(jvertex) != _global_boundary.end() )
      //     return true;
      return false;
    }
  }

}

}  // end namespace mesh
