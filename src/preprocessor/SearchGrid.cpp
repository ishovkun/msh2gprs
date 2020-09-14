#include "SearchGrid.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include <queue>

namespace gprs_data {

using Location = std::array<int, 3>;

SearchGrid::SearchGrid(const mesh::Mesh & grid)
    : _grid(grid)
{
  compute_stepping_and_find_origin_();
  _mapping.resize(n_cells());
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    map_cell_(cell->index());
}

void SearchGrid::map_cell_(size_t cell_index)
{
 std::queue<size_t> search_queue;
 const auto cell_poly =_grid.cell(cell_index).polyhedron();
 const auto c = cell_poly->center();
 search_queue.push(global_index(index(c)));
 angem::CollisionGJK<double> collision;
 std::set<size_t> processed;
 while (!search_queue.empty())
 {
   const size_t search_cell = search_queue.front();
   search_queue.pop();
   const auto search_poly = get_cell_hex_(search_cell);
   if (collision.check(*cell_poly, search_poly))
   {
     _mapping[search_cell].push_back(cell_index);
     for (size_t neighbor : neighbors(search_cell))
       if (processed.find(neighbor) == processed.end())
       {
         processed.insert(neighbor);
         search_queue.push(neighbor);
       }
   }
 }
}

size_t SearchGrid::global_index(const std::array<int, 3> &loc) const
{
  return loc[2] * _dims[0] * _dims[1] + loc[1] * _dims[0] + loc[0];
}

size_t SearchGrid::global_index(int i, int j, int k) const
{
  return k*_dims[0]*_dims[1] + j*_dims[0] + i;
}

std::array<int, 3> SearchGrid::local_index(size_t idx) const
{
  if (idx >= n_cells()) throw std::invalid_argument(std::to_string(idx) + " out of bounds");
  size_t a = idx % (_dims[0]*_dims[1]);
  Location loc;
  loc[0] = a % _dims[0];
  loc[1] = (a-loc[0]) / _dims[0];
  loc[2] = (idx - loc[0] - loc[1]*_dims[0]) / _dims[0] / _dims[1];
  return loc;
}


angem::Hexahedron<double> SearchGrid::get_cell_hex_(size_t idx) const
{
  const auto loc = local_index(idx);
  std::vector<angem::Point<3,double>> verts;
  verts.emplace_back(_origin[0] + loc[0]*_stepping[0],
                     _origin[1] + loc[1]*_stepping[1],
                     _origin[2] + loc[2]*_stepping[2]);  // v0

  verts.emplace_back(_origin[0] + (loc[0]+1)*_stepping[0],
                     _origin[1] + loc[1]*_stepping[1],
                     _origin[2] + loc[2]*_stepping[2]);  // v1

  verts.emplace_back(_origin[0] + (loc[0]+1)*_stepping[0],
                     _origin[1] + (loc[1]+1)*_stepping[1],
                     _origin[2] + loc[2]*_stepping[2]);  // v2

  verts.emplace_back(_origin[0] + loc[0]*_stepping[0],
                     _origin[1] + (loc[1]+1)*_stepping[1],
                     _origin[2] + loc[2]*_stepping[2]);  // v3

  verts.emplace_back(_origin[0] + loc[0]*_stepping[0],
                     _origin[1] + loc[1]*_stepping[1],
                     _origin[2] + (loc[2]+1)*_stepping[2]);  // v4

  verts.emplace_back(_origin[0] + (loc[0]+1)*_stepping[0],
                     _origin[1] + loc[1]*_stepping[1],
                     _origin[2] + (loc[2]+1)*_stepping[2]);  // v5

  verts.emplace_back(_origin[0] + (loc[0]+1)*_stepping[0],
                     _origin[1] + (loc[1]+1)*_stepping[1],
                     _origin[2] + (loc[2]+1)*_stepping[2]);  // v6

  verts.emplace_back(_origin[0] + loc[0]*_stepping[0],
                     _origin[1] + (loc[1]+1)*_stepping[1],
                     _origin[2] + (loc[2]+1)*_stepping[2]);  // v7
  std::vector<size_t> indices(8);
  std::iota(indices.begin(), indices.end(), 0);
  return angem::Hexahedron<double>(verts, indices);
}

bool SearchGrid::in_bounds(const angem::Point<3,double> & p) const noexcept
{
  for (size_t i = 0; i < 3; ++i)
    if (p[i] < _origin[i] || p[i] > _origin[i] + _stepping[i] * _dims[i])
      return false;
  return true;
}

std::array<int, 3> SearchGrid::index(const angem::Point<3,double> & p) const
{
  if (!in_bounds(p))
    throw std::invalid_argument("point " + std::to_string(p[0]) + " " +
                                std::to_string(p[1]) + " " +
                                std::to_string(p[2]) + " not in bounds");
  Location loc;
  for (size_t i = 0; i < 3; ++i)
    loc[i] = (size_t)((p[i] - _origin[i]) / _stepping[i]);
  return loc;
}


size_t SearchGrid::n_cells() const noexcept
{
  return std::accumulate(_dims.begin(), _dims.end(), 1, std::multiplies<size_t>());
}

void SearchGrid::compute_stepping_and_find_origin_()
{
  // origin is the point with min(x), min(y), min(z) coordinates
  // stepping is the average grid spacing in x, y , and z directions
  _origin[0] = _origin[1] = _origin[2] = std::numeric_limits<double>::max();
  angem::Point<3,double> corner;
  corner[0] = corner[1] = corner[2] =  std::numeric_limits<double>::lowest();
  assert( _origin[0] == std::numeric_limits<double>::max());
  _stepping = {0,0,0};
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    angem::Point<3,double> min_point, max_point;
    for (const size_t vertex : cell->vertices())
    {
      const auto & v = _grid.vertex(vertex);
      for (size_t i = 0; i < 3; ++i)
      {
        min_point[i] = std::min(min_point[i], v[i]);
        max_point[i] = std::max(max_point[i], v[i]);
        _origin[i] = std::min(_origin[i], min_point[i]);
        corner[i] = std::max(corner[i], max_point[i]);
      }
      _stepping += max_point - min_point;
    }
  }
  _stepping /= _grid.n_active_cells();
  for (size_t i = 0; i < 3; ++i)
  {
    _dims[i] = (size_t)((corner[i] - _origin[i]) / _stepping[i]);
    _stepping[i] = (corner[i] - _origin[i]) / _dims[i];
  }
}

std::vector<size_t> SearchGrid::neighbors(size_t search_cell) const
{
  const auto loc = local_index(search_cell);
  std::vector<size_t> result;
  result.reserve(6);
  add_neighbor_(loc[0]+1, loc[1], loc[2], result);
  add_neighbor_(loc[0]-1, loc[1], loc[2], result);
  add_neighbor_(loc[0]  , loc[1]+1, loc[2], result);
  add_neighbor_(loc[0]  , loc[1]-1, loc[2], result);
  add_neighbor_(loc[0]  , loc[1]  , loc[2]+1, result);
  add_neighbor_(loc[0]  , loc[1]  , loc[2]-1, result);
  return result;
}

void SearchGrid::add_neighbor_(int i, int j, int k, std::vector<size_t> & dst) const
{
  if (in_bounds(i, j, k))
    dst.push_back(global_index(i,j,k));
}


bool SearchGrid::in_bounds(int i, int j, int k) const noexcept
{
  if (i > 0 && j > 0 && k > 0 &&
      i < _dims[0] && j < _dims[1] && k < _dims[2])
    return true;
  else return false;
}


}  // end namespace gprs_data
