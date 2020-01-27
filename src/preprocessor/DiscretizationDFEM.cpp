#include "DiscretizationDFEM.hpp"
#include <numeric>  // provides std;:iota

namespace gprs_data
{
using Point = angem::Point<3,double>;

DiscretizationDFEM::DiscretizationDFEM(const mesh::Mesh & grid)
    : _grid(grid)
{}

#ifdef WITH_GMSH

void DiscretizationDFEM::build()
{
  // build_gmsh_simple();
  // exit(0);
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    std::cout << "cell.index() = " << cell->index() << std::endl;
    const auto poly = cell->polyhedron();
    build_(*poly);
    exit(0);
  }
}

void DiscretizationDFEM::build_(const angem::Polyhedron<double> & cell)
{
  gmsh::initialize();
  build_grid_(cell);
  build_shape_functions_();
  gmsh::finalize();
}

void DiscretizationDFEM::build_grid_(const angem::Polyhedron<double> & cell) const
{
  gmsh::option::setNumber("General.Terminal", 1);

  const double discr_element_size = 0.1 * compute_element_size_(cell);
  std::cout << "discr_element_size = " << discr_element_size << std::endl;
  // build points
  const std::vector<Point> & vertices = cell.get_points();
  std::cout << "n_vertices = " << vertices.size() << std::endl;
  for (std::size_t i=0; i < vertices.size(); ++i)
  {
    const Point & vertex = vertices[i];
    gmsh::model::geo::addPoint(vertex.x(), vertex.y(), vertex.z(),
                               discr_element_size, /*tag = */ i);
  }

  // build lines (edges)
  const auto edges = cell.get_edges();
  for (std::size_t i=0; i<edges.size(); ++i)
  {
    const std::pair<size_t,size_t> & edge = edges[i];
    gmsh::model::geo::addLine(edge.first, edge.second, i);
    // gmsh::model::addPhysicalGroup(1, {i}, i);
  }

  // build faces
  const std::vector<std::vector<std::size_t>> & faces = cell.get_faces();
  std::vector<angem::Polygon<double>> face_polys = cell.get_face_polygons();
  for (std::size_t i=0; i<face_polys.size(); ++i)
  {
    const auto & face = faces[i];
    const auto & poly = face_polys[i];
    std::vector<int> edge_markers;
    for (const angem::Edge & face_edge : poly.get_edges())
    {
      const std::pair<size_t,size_t> cell_edge_ordered = {face[face_edge.first],
                                                         face[face_edge.second]};
      const std::pair<size_t,size_t> cell_edge_unordered = std::minmax(face[face_edge.first],
                                                                       face[face_edge.second]);
      // get edge marker
      const auto it_edge = std::find_if( edges.begin(), edges.end(),
                    [&cell_edge_unordered](const auto & it)->bool
                    {
                      return it.first == cell_edge_unordered.first &&
                          it.second == cell_edge_unordered.second;
                    });
      assert(it_edge != edges.end());;
      const int edge_marker = static_cast<int>(std::distance(edges.begin(), it_edge));
      // figure out the sign of the tage
      if (cell_edge_ordered.first == cell_edge_unordered.first)
        edge_markers.push_back(edge_marker);
      else  // inverse orientation
        edge_markers.push_back(-edge_marker);
    }
    // create line loop and surface
    // NOTE: curve and surface loop must start from 1, otherwise gmsh
    // throws an error, ergo i+1
    gmsh::model::geo::addCurveLoop(edge_markers, static_cast<int>(i+1));
    gmsh::model::geo::addPlaneSurface({static_cast<int>(i+1)}, static_cast<int>(i+1));
  }

  // gmsh::model::geo::synchronize();
  // gmsh::model::mesh::generate(2);
  // gmsh::write("cell.msh");
  // gmsh::finalize();

  // create volume from surfaces
  std::vector<int> surfaces(faces.size());
  std::iota(surfaces.begin(), surfaces.end(), 1);
  gmsh::model::geo::addSurfaceLoop(surfaces, 1);
  gmsh::model::geo::addVolume({1}, 1);

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);
  gmsh::write("cell.msh");
}

double DiscretizationDFEM::compute_element_size_(const angem::Polyhedron<double> & cell) const
{
  double min_edge_length = std::numeric_limits<double>::max();
  const std::vector<Point> & vertices = cell.get_points();
  for (const auto & edge : cell.get_edges())
  {
    const double h = vertices[edge.first].distance(vertices[edge.second]);
    min_edge_length = std::min( min_edge_length, h);
  }
  return min_edge_length;
}

void DiscretizationDFEM::build_shape_functions_()
{
  std::vector<int> element_types;
  std::vector<std::vector<std::size_t> > element_tags;
  std::vector<std::vector<std::size_t> > node_tags;
  gmsh::model::mesh::getElements(element_types, element_tags, node_tags);
  std::cout << "tags" << std::endl;
  for (auto type : element_types)
    std::cout << type << std::endl;

  // 1. make sure that we only have tetras (look at element_types)
  // 2. ask gmsh to provide gaussian points, shape functions, and jacobians
  // for our tetras
  // 3. build element jacobian for the homogeneous laplace equation
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace gprs_data
