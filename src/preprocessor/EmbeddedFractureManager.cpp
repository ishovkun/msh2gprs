#include "EmbeddedFractureManager.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include "angem/Collisions.hpp"    // angem::split

namespace gprs_data {

using std::vector;
using Point = angem::Point<3,double>;
const int MARKER_BELOW_FRAC = 0;
const int MARKER_ABOVE_FRAC = 1;
const int MARKER_FRAC = 2;

EmbeddedFractureManager::
EmbeddedFractureManager(const std::vector<EmbeddedFractureConfig> &config,
                        const EDFMMethod edfm_method,
                        SimData & data)
    : config(config), m_method(edfm_method), data(data), m_split_grid(data.grid)
{}

void EmbeddedFractureManager::split_cells()
{
  for (auto & frac : config)  // non-const since we can shift it
  {
    vector<size_t> cells;
    // iteratively shift fracture if it collides with any grid vertices
    size_t iter = 0;
    while (!find_edfm_cells_(*frac.body, cells) && iter < 100)
    {
      cells.clear();
      iter++;
    }

    split_cells_(*frac.body, cells);
  }
}

void EmbeddedFractureManager::split_cells_(angem::Polygon<double> & fracture,
                                           std::vector<size_t> & cells)
{
  const auto & plane = fracture.plane;
  for (const size_t icell : cells)
  {
    std::cout << "icell = " << icell << std::endl;
    mesh::Cell & old_cell = m_split_grid.cell(icell);
    const std::unique_ptr<angem::Polyhedron<double>> polyhedron = old_cell.polyhedron();

    // Bookkeeping:
    //  fill polygroup's internal set with the existing vertex coordinates
    // in order to have a map of those to the global vertex indices,
    // which will come in handy when inserting new splitted cells into grid.
    // We can do it because splitting will insert the same vertices plus
    // those that appeared due to plase-face intersection.
    angem::PolyGroup<double> split;
    std::vector<size_t> global_vertex_indices;
    for (const Point & p : polyhedron->get_points())
      global_vertex_indices.push_back(split.vertices.insert(p));

    std::cout << "old_size = " << global_vertex_indices.size() << std::endl;
    angem::split(*polyhedron, plane, split, MARKER_BELOW_FRAC,
                 MARKER_ABOVE_FRAC, MARKER_FRAC);
    std::cout << "new_size = " << split.vertices.size() << std::endl;
    // check we actually split something
    assert( split.vertices.size() > global_vertex_indices.size() );

    // insert new vertices
    for (size_t i = global_vertex_indices.size(); i < split.vertices.size(); ++i)
    {
      const size_t new_vertex_index = m_split_grid.n_vertices();
      m_split_grid.vertices().push_back(split.vertices[i]);
      global_vertex_indices.push_back(new_vertex_index);
    }

    // make two polyhedra that will form the new vertices
    vector<size_t> polyhedra_above_faces, polyhedra_below_faces;
    for (size_t i = 0; i < split.polygons.size(); i++)
    {
      if ( split.markers[i] == MARKER_BELOW_FRAC || split.markers[i] == MARKER_FRAC )
        polyhedra_below_faces.push_back(split.polygons[i]);
      if ( split.markers[i] == MARKER_ABOVE_FRAC || split.markers[i] == MARKER_FRAC )
        polyhedra_above_faces.push_back(split.polygons[i]);
    }

    exit(0);
  }
}

bool EmbeddedFractureManager::find_edfm_cells_(angem::Polygon<double> & fracture,
                                               std::vector<size_t> & cells) const
{
  const auto & grid = data.grid;
  // performs fast collision check
  angem::CollisionGJK<double> collision;

  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    const std::unique_ptr<angem::Polyhedron<double>> polyhedron = cell->polyhedron();
    // const angem::Polyhedron<double> & poly_cell = *pol;
    if (collision.check(fracture, *polyhedron))
    {
      cells.push_back(cell->index());

      // check if some vertices are too close to the fracture
      // and if so move a fracture a little bit
      const std::vector<Point> & vertices = polyhedron->get_points();
      for (const Point & vertex : vertices)
      {
        const double dist_vert_center = (polyhedron->center() - vertex).norm();
        if ( std::fabs( fracture.plane.distance(vertex) / dist_vert_center) < 1e-4 )
        {
          const double h = vertices[1].distance(vertices[0]);
          const Point shift = h/5 * fracture.plane.normal();
          fracture.move(shift);
          return false;
        }
        
      }

    }
  }

  return true;
}

}  // end namespace gprs_data
