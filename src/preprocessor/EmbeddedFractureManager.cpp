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
    m_split_grid.split_cell(old_cell, plane);
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
        if ( std::fabs( fracture.plane.signed_distance(vertex) / dist_vert_center) < 1e-4 )
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
