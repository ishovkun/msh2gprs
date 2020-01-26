#include "DiscretizationFEM.hpp"
#include <numeric>  // provides std;:iota

namespace gprs_data
{
using Point = angem::Point<3,double>;

DiscretizationFEM::DiscretizationFEM(const mesh::Mesh & grid)
    : _grid(grid)
{}

#ifdef WITH_GMSH

void build_gmsh_simple()
{
  // Before using any functions in the C++ API, Gmsh must be initialized.
  gmsh::initialize();

  // By default Gmsh will not print out any messages: in order to output
  // messages on the terminal, just set the standard Gmsh option
  // "General.Terminal" (same format and meaning as in .geo files) using
  // gmsh::option::setNumber():
  gmsh::option::setNumber("General.Terminal", 1);

  // This adds a new model, named "t1". If gmsh::model::add() is not called, a
  // new default (unnamed) model will be created on the fly, if necessary.
  gmsh::model::add("t1");

  // The C++ API provides direct access to the internal CAD kernels. The
  // built-in CAD kernel was used in t1.geo: the corresponding API functions
  // live in the "gmsh::model::geo" namespace. To create geometrical points with
  // the built-in CAD kernel, one thus uses gmsh::model::geo::addPoint():
  //
  // - the first 3 arguments are the point coordinates (x, y, z)
  //
  // - the next (optional) argument is the target mesh size close to the point
  //
  // - the last (optional) argument is the point tag
  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
  gmsh::model::geo::addPoint(.1, 0,  0, lc, 2);
  gmsh::model::geo::addPoint(.1, .3, 0, lc, 3);
  gmsh::model::geo::addPoint(0,  .3, 0, lc, 4);

  // The API to create lines with the built-in kernel follows the same
  // conventions: the first 2 arguments are point tags, the last (optional one)
  // is the line tag.
  gmsh::model::geo::addLine(1, 2, 1);
  gmsh::model::geo::addLine(3, 2, 2);
  gmsh::model::geo::addLine(3, 4, 3);
  gmsh::model::geo::addLine(4, 1, 4);

  // The philosophy to construct curve loops and surfaces is similar: the first
  // argument is now a vector of integers.
  gmsh::model::geo::addCurveLoop({4, 1, -2, 3}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);

  // Physical groups are defined by providing the dimension of the group (0 for
  // physical points, 1 for physical curves, 2 for physical surfaces and 3 for
  // phsyical volumes) followed by a vector of entity tags. The last (optional)
  // argument is the tag of the new group to create.
  gmsh::model::addPhysicalGroup(0, {1, 2}, 1);
  gmsh::model::addPhysicalGroup(1, {1, 2}, 2);
  gmsh::model::addPhysicalGroup(2, {1}, 6);

  // Physical names are also defined by providing the dimension and tag of the
  // entity.
  gmsh::model::setPhysicalName(2, 6, "My surface");

  // Before it can be meshed, the internal CAD representation must be
  // synchronized with the Gmsh model, which will create the relevant Gmsh data
  // structures. This is achieved by the gmsh::model::geo::synchronize() API
  // call for the built-in CAD kernel. Synchronizations can be called at any
  // time, but they involve a non trivial amount of processing; so while you
  // could synchronize the internal CAD data after every CAD command, it is
  // usually better to minimize the number of synchronization points.
  gmsh::model::geo::synchronize();

  // We can then generate a 2D mesh...
  gmsh::model::mesh::generate(2);

  // ... and save it to disk
  gmsh::write("t1.msh");

  // Remember that by default, if physical groups are defined, Gmsh will export
  // in the output mesh file only those elements that belong to at least one
  // physical group. To force Gmsh to save all elements, you can use
  //
  gmsh::option::setNumber("Mesh.SaveAll", 1);

  // This should be called at the end:
  gmsh::finalize();
}

void DiscretizationFEM::build()
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

void DiscretizationFEM::build_(const angem::Polyhedron<double> & cell)
{
  gmsh::initialize();
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
    std::cout << "face " << i << ": ";
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
    // std::cout << cell_edge_ordered.first << "-"<<cell_edge_ordered.second << " ";
    }
    for (auto m : edge_markers)
      std::cout << m << " ";
    std::cout << std::endl;
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
  gmsh::finalize();
}

double DiscretizationFEM::compute_element_size_(const angem::Polyhedron<double> & cell)
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

#else
void DiscretizationFEM::build() {}
#endif

}  // end namepsace gprs_data
