#ifdef WITH_EIGEN
#include "TributaryRegion3dVertices.hpp"
#include "VTKWriter.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

std::vector<mesh::Edge> get_edges_with_vertex(const size_t vertex,
                                              const mesh::Face & face)
{
  const auto & vertices = face.vertices();
  const size_t ivertex = std::distance(vertices.begin(), std::find(
      vertices.begin(), vertices.end(), vertex));
  size_t iv1, iv2;
  if (ivertex < vertices.size() - 1)
      iv2 = ivertex + 1;
  else iv2 = 0;
  if (ivertex > 0)
      iv1 = ivertex - 1;
  else iv1 = vertices.size() - 1;

  return {{vertex, vertices[iv1]}, {vertex, vertices[iv2]}};
}

TributaryRegion3dVertices::TributaryRegion3dVertices(PolyhedralElementBase & element)
    : TributaryRegion3dBase(element)
{
  const auto & cell = _element._parent_cell;
  const auto & vertices = cell.vertices();

  for (size_t ivertex = 0; ivertex < vertices.size(); ++ivertex)
  {
    std::vector<std::vector<size_t>> tributary_cell_faces;
    angem::PointSet<3, double> tributary_cell_vertices;
    const size_t vertex = vertices[ivertex];
    for (const auto face : cell.faces())
    {
      if (face->has_vertex(vertex))
      {
        const auto edges = get_edges_with_vertex(vertex, *face);
        build_tributary_cell_faces_(edges, *face, tributary_cell_faces, tributary_cell_vertices);
      }
    }
    _tributary.emplace_back(tributary_cell_vertices.points, tributary_cell_faces);
    IO::VTKWriter::write(_tributary.back(), "tributary-" + std::to_string(cell.index()) +
                     "-" + std::to_string(ivertex) + ".vtk");
  }
  exit(0);
}

void TributaryRegion3dVertices::
build_tributary_cell_faces_(const std::vector<mesh::Edge> & edges,
                            const mesh::Face & face,
                            std::vector<std::vector<size_t>> & tributary_faces,
                            angem::PointSet<3,double> & tributary_vertices)
{
  const auto & grid = _element._parent_grid;
  const size_t v0 = tributary_vertices.insert(grid.vertex(edges[0].first));
  const Point p1 = 0.5 * (grid.vertex(edges[0].first)  + grid.vertex(edges[0].second));
  const Point p2 = 0.5 * (grid.vertex(edges[1].first)  + grid.vertex(edges[1].second));
  const size_t v1 = tributary_vertices.insert(p1);
  const size_t v2 = tributary_vertices.insert(p2);
  const size_t v3 = tributary_vertices.insert(face.center());
  tributary_faces.push_back({v0, v1, v3, v2});
  // cell center
  const size_t vc = tributary_vertices.insert(_element._parent_cell.center());
  tributary_faces.push_back({v0, v1, v2});
  // tributary_faces.push_back({vc, v1, v2});
}

}  // end namespace discretization

#endif
