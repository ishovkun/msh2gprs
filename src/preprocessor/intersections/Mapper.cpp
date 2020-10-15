#include "Mapper.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include <queue>
#include "logger/Logger.hpp"

namespace gprs_data {

using Location = std::array<int, 3>;

Mapper::Mapper(const mesh::Mesh & grid)
    : _grid(grid), _cartesian(build_grid_()), _num_cells_cached(_grid.n_cells_total())
{
  _cartesian.log_stats();
  _mapping.resize(_cartesian.n_cells());
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    map_cell(cell->index());
}

void Mapper::map_cell(size_t cell_index)
{
 std::queue<size_t> search_queue;
 const auto cell_poly =_grid.cell(cell_index).polyhedron();
 // start the cearch in the voxel that contains grid cell center
 search_queue.push(_cartesian.find_cell(cell_poly->center()));
 angem::CollisionGJK<double> collision;
 std::set<size_t> touched;
 while (!search_queue.empty())
 {
   const size_t search_cell = search_queue.front();
   search_queue.pop();
   const auto search_poly = _cartesian.get_voxel(search_cell);
   if (collision.check(*cell_poly, search_poly))
   {
     _mapping[search_cell].push_back(cell_index);
     for (size_t const neighbor : _cartesian.neighbors(search_cell))
       if (touched.find(neighbor) == touched.end())
       {
         touched.insert(neighbor);
         search_queue.push(neighbor);
       }
   }
 }
}

UniformCartesianGrid Mapper::build_grid_() const
{
  // origin is the point with min(x), min(y), min(z) coordinates
  // stepping is the average grid spacing in x, y , and z directions
  const double max_double = std::numeric_limits<double>::max();
  const double min_double = std::numeric_limits<double>::lowest();
  angem::Point<3,double> loc_min, loc_max, glob_min, glob_max, stepping;
  glob_min[0] = glob_min[1] = glob_min[2] = max_double;
  glob_max[0] = glob_max[1] = glob_max[2] = min_double;
  stepping = {0,0,0};
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    loc_min[0] = loc_min[1] = loc_min[2] = max_double;
    loc_max[0] = loc_max[1] = loc_max[2] = min_double;
    for (const size_t vertex : cell->vertices())
    {
      const auto & v = _grid.vertex(vertex);
      for (size_t i = 0; i < 3; ++i)
      {
        loc_min[i] = std::min(loc_min[i], v[i]);
        loc_max[i] = std::max(loc_max[i], v[i]);
      }
    }
    stepping += loc_max - loc_min;
    for (size_t i = 0; i < 3; ++i)
    {
      glob_min[i] = std::min(glob_min[i], loc_min[i]);
      glob_max[i] = std::max(glob_max[i], loc_max[i]);
    }
  }
  stepping /= _grid.n_active_cells();
  std::array<size_t,3> dims;
  for (size_t i = 0; i < 3; ++i)
  {
    dims[i] = (size_t)((glob_max[i] - glob_min[i]) /(double) stepping[i]);
    stepping[i] = (glob_max[i] - glob_min[i]) / dims[i];
  }
  return UniformCartesianGrid(glob_min, stepping, dims);
}

const std::list<size_t> & Mapper::mapping(size_t search_cell)
{
  if (_num_cells_cached < _grid.n_cells_total())
    update_mapping();
  if (search_cell < _mapping.size())
    return _mapping[search_cell];
  else throw std::invalid_argument("invalid cell " + std::to_string(search_cell));
}

void Mapper::update_mapping()
{
  for (size_t icell = _num_cells_cached; icell < _grid.n_cells_total(); ++icell)
    map_cell(icell);
  _num_cells_cached = _grid.n_cells_total();
}


}  // end namespace gprs_data
