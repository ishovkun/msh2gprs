#ifdef WITH_EIGEN
#include "TributaryRegion3dVertices.hpp"
#include "VTKWriter.hpp"      // provides IO::VTKWriter
#include "../EdgeComparison.hpp"  // provides Edge and EdgeComparison


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
  const auto & grid = _element._parent_grid;
  const auto & cell = _element._parent_cell;
  const auto & vertices = cell.vertices();

  const auto parent_vertex_faces = _element.map_parent_vertices_to_parent_faces_();
  const auto pair_markers_to_edge = edgecmp::EdgeComparison::get_edges( parent_vertex_faces,
                                                                        _element._parent_cell.polyhedron()->get_edges());
  const auto parent_faces = cell.faces();
  const Point c = cell.center();
  // std::cout << "parent_faces.size() = " << parent_faces.size() << std::endl;
  // IO::VTKWriter::write_geometry(grid, cell, "parent.vtk");
  

  for (size_t ivertex = 0; ivertex < vertices.size(); ++ivertex)
  {
    // these two entities store information for the faces of a tributary volume
    angem::PointSet<3, double> tributary_cell_vertices;
    std::vector<std::vector<size_t>> tributary_cell_faces;
    const size_t v = tributary_cell_vertices.insert(grid.vertex(vertices[ivertex]));
    const size_t cc = tributary_cell_vertices.insert(c);
    // std::cout << "ivertex = " << ivertex << std::endl;

    // loop pairs of faces adjacent to the vertex
    const auto & markers = parent_vertex_faces[ivertex];
    for (auto it1 = markers.begin(); it1 != markers.end(); ++it1)
      for (auto it2 = std::next(it1, 1); it2 != markers.end(); ++it2)
        if (pair_markers_to_edge.contains(*it1, *it2))
        {
          // std::cout << "\t" << *it1 << " " << *it2 << std::endl;
          const auto * face1 = parent_faces[*it1 - 1];
          const auto * face2 = parent_faces[*it2 - 1];
          const size_t f1c = tributary_cell_vertices.insert(face1->center());
          const size_t f2c = tributary_cell_vertices.insert(face2->center());

          // there are two edges that separate two faces adjacent to ivertex
          // in the case of a handing node; otherwise 1 edge
          const auto & edges = pair_markers_to_edge.get_data(*it1, *it2);
          // if (edges.size() == 1)  // not a hanging node
          for (const auto & edge : edges)
          {
            // const auto & edge = edges[0];
            const size_t jvertex = (edge.either() == ivertex) ? edge.other(edge.either()) :
                                                                edge.either();
            const Point edge_center = 0.5 * (grid.vertex(vertices[ivertex]) +
                                             grid.vertex(vertices[jvertex]));
            const size_t ec = tributary_cell_vertices.insert(edge_center);
            tributary_cell_faces.push_back({v, ec, f1c});
            tributary_cell_faces.push_back({ec, f1c, cc});
            tributary_cell_faces.push_back({v, ec, f2c});
            tributary_cell_faces.push_back({ec, f2c, cc});
          }
          // else if (vertices.size() == 2)  // not a hanging node
          // {
          //   throw std::runtime_error("hanging not not implemented yet");
          // }
          // else
          //   throw std::runtime_error("unexpected situation; should not happen");
        }
    _tributary.emplace_back(tributary_cell_vertices.points, tributary_cell_faces);
    // IO::VTKWriter::write(_tributary.back(), "tributary-" + std::to_string(cell.index()) +
    //                  "-" + std::to_string(ivertex) + ".vtk");
  }
  // exit(0);
}

void TributaryRegion3dVertices::
build_tributary_cell_faces_(const std::vector<mesh::Edge> & edges,
                            const mesh::Face & face,
                            std::vector<std::vector<size_t>> & tributary_faces,
                            angem::PointSet<3,double> & tributary_vertices)
{
  const auto & grid = _element._parent_grid;
  const size_t v0 = tributary_vertices.insert(grid.vertex(edges[0].first));
  const Point p1 = 0.5 * (grid.vertex(edges[0].first) + grid.vertex(edges[0].second));
  const Point p2 = 0.5 * (grid.vertex(edges[1].first) + grid.vertex(edges[1].second));
  const size_t v1 = tributary_vertices.insert(p1);
  const size_t v2 = tributary_vertices.insert(p2);
  const size_t v3 = tributary_vertices.insert(face.center());
  tributary_faces.push_back({v0, v1, v3, v2});
  // cell center
  const size_t vc = tributary_vertices.insert(_element._parent_cell.center());
  tributary_faces.push_back({v0, v1, v2});
  // tributary_faces.push_back({vc, v1, v2});
}

void TributaryRegion3dVertices::
build_face_(const size_t vert1, const size_t vert2,
            const mesh::Face & face,
            std::vector<std::vector<size_t>> & tributary_faces,
            angem::PointSet<3,double> & tributary_vertices) const
{

}

}  // end namespace discretization

#endif
