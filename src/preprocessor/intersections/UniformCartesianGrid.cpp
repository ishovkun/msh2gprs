#include "UniformCartesianGrid.hpp"
#include "logger/Logger.hpp"
#include <algorithm>  // clamp

namespace gprs_data {

UniformCartesianGrid::UniformCartesianGrid(const angem::Point<3,double> &origin,
                                           const angem::Point<3,double> &stepping,
                                           const std::array<size_t,3> &dims)
    : _origin(origin), _stepping(stepping), _dims(dims)
{
  if (std::any_of(_dims.begin(), _dims.end(), [](size_t value) { return value == 0; }))
    throw std::invalid_argument("Wrong grid dimensions");
}

void UniformCartesianGrid::log_stats() const noexcept
{
  logging::debug() << "Indexing, overlay grid: " << n_cells()
                   << " dimensions = (" << _dims[0] << ", " << _dims[1] << ", "
                   << _dims[2] << ")" << std::endl;
  logging::debug() << "bounding box " << _origin << " to "
                   << _origin[0] + _stepping[0]*_dims[0] << " "
                   << _origin[1] + _stepping[1]*_dims[1] << " "
                   << _origin[2] + _stepping[2]*_dims[2] << std::endl;
}

size_t UniformCartesianGrid::n_cells() const noexcept
{
  return _dims[0] * _dims[1] * _dims[2];
}

size_t UniformCartesianGrid::cell_index(int i, int j, int k) const
{
  return k*_dims[0]*_dims[1] + j*_dims[0] + i;
}

std::array<int,3> UniformCartesianGrid::get_ijk(size_t idx) const
{
  const size_t rem = idx % (_dims[0] * _dims[1]);
  const int i = rem % _dims[0];
  const int j = (rem - i) / _dims[0];
  const int k = (idx - i - _dims[0]*j) / _dims[0] / _dims[1];
  return {i, j, k};
}

bool UniformCartesianGrid::is_valid_cell(int i, int j, int k) const noexcept
{
  if (i >= 0 && j >= 0 && k >= 0 && i < _dims[0] && j < _dims[1] && k < _dims[2])
    return true;
  else return false;
}

size_t UniformCartesianGrid::find_cell(const angem::Point<3,double> & p) const
{
  if (!in_bounds(p))
    throw std::invalid_argument("point " + std::to_string(p[0]) + " " +
                                std::to_string(p[1]) + " " +
                                std::to_string(p[2]) + " not in bounds");
  std::array<int,3> loc;
  for (size_t i = 0; i < 3; ++i)
  {
    loc[i] = (int) std::floor((p[i] - _origin[i]) / _stepping[i]);
    loc[i] = std::clamp(loc[i], int(0), int(_dims[i]) - 1);
  }

  return cell_index(loc[0], loc[1], loc[2]);
}

bool UniformCartesianGrid::in_bounds(const angem::Point<3,double> & p) const noexcept
{
  for (size_t i = 0; i < 3; ++i)
    if (p[i] < _origin[i] || p[i] > _origin[i] + _stepping[i] * _dims[i])
      return false;
  return true;
}

angem::Hexahedron<double> UniformCartesianGrid::get_voxel(size_t idx) const
{
  if (idx >= n_cells()) throw std::invalid_argument("cell does not exist");
  const auto ijk = get_ijk(idx);
  const int i = ijk[0];
  const int j = ijk[1];
  const int k = ijk[2];
  const auto v0 = vertex(i,   j,   k);
  const auto v1 = vertex(i+1, j,   k);
  const auto v2 = vertex(i+1, j+1, k);
  const auto v3 = vertex(i,   j+1, k);
  const auto v4 = vertex(i,   j,   k+1);
  const auto v5 = vertex(i+1, j,   k+1);
  const auto v6 = vertex(i+1, j+1, k+1);
  const auto v7 = vertex(i,   j+1, k+1);
  const std::vector<angem::Point<3,double>> vertices = {v0, v1, v2, v3, v4, v5, v6, v7};
  std::vector<size_t> indices(8);
  std::iota(indices.begin(), indices.end(), 0);
  return angem::Hexahedron<double>(vertices,indices);
}

angem::Point<3,double> UniformCartesianGrid::vertex(int i, int j, int k)  const
{
  if (is_valid_vertex(i, j, k))
  {
    angem::Point<3,double> result;
    std::array<int, 3> ijk{i, j, k};
    for (size_t d = 0; d < 3; ++d)
      result[d] = _origin[d] + _stepping[d] * ijk[d];
    return result;
  }
  else throw std::invalid_argument("Vertex is invalid " + std::to_string(i) + " " +
                                   std::to_string(j) + " " + std::to_string(k));
}

bool UniformCartesianGrid::is_valid_vertex(int i, int j, int k) const noexcept
{
  if (i >= 0 && j >= 0 && k >= 0 && i < _dims[0]+1 && j < _dims[1]+1 && k < _dims[2]+1)
    return true;
  else return false;
}

std::vector<size_t> UniformCartesianGrid::neighbors(size_t cell_idx) const
{
  std::vector<size_t> result;
  result.reserve(6);
  const auto loc = get_ijk(cell_idx);
  add_neighbor_(loc[0]+1, loc[1], loc[2], result);
  add_neighbor_(loc[0]-1, loc[1], loc[2], result);
  add_neighbor_(loc[0]  , loc[1]+1, loc[2], result);
  add_neighbor_(loc[0]  , loc[1]-1, loc[2], result);
  add_neighbor_(loc[0]  , loc[1]  , loc[2]+1, result);
  add_neighbor_(loc[0]  , loc[1]  , loc[2]-1, result);
  return result;
}

void UniformCartesianGrid::add_neighbor_(int i, int j, int k, std::vector<size_t> & dst) const
{
  if (is_valid_cell(i, j, k))
    dst.push_back(cell_index(i,j,k));
}

angem::Hexahedron<double> UniformCartesianGrid::get_bounding_box() const noexcept
{
  const size_t nvx = _dims[0] + 1;
  const size_t nvy = _dims[1] + 1;
  const size_t nvz = _dims[2] + 1;
  const auto v0 = vertex(0,       0,       0);
  const auto v1 = vertex(nvx - 1, 0,       0);
  const auto v2 = vertex(nvx - 1, nvy - 1, 0);
  const auto v3 = vertex(0,       nvy - 1, 0);
  const auto v4 = vertex(0,       0,       nvz - 1);
  const auto v5 = vertex(nvx - 1, 0,       nvz - 1);
  const auto v6 = vertex(nvx - 1, nvy - 1, nvz - 1);
  const auto v7 = vertex(0,       nvy - 1, nvz - 1);
  const std::vector<angem::Point<3,double>> vertices = {v0, v1, v2, v3, v4, v5, v6, v7};
  std::vector<size_t> indices(8);
  std::iota(indices.begin(), indices.end(), 0);
  return angem::Hexahedron<double>(vertices,indices);
}

double UniformCartesianGrid::stepping(const size_t dir) const
{
  if (dir >= 3) throw std::invalid_argument("wrong direction (only 0,1,2 are allowed)");
  return _stepping[dir];
}


}  // end namespace gprs_data
