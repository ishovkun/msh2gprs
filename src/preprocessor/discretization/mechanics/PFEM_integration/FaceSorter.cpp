#include "FaceSorter.hpp"

namespace discretization {

FaceSorter::FaceSorter(const mesh::Face & parent_face, const mesh::Mesh & subgrid) noexcept
    : _grid(subgrid),
      _pmap(point_comparator(angem::Plane<double>(parent_face.vertex_coordinates())))
{}

void FaceSorter::sort(std::vector<size_t> & face_indices)
{
  for (const size_t iface : face_indices)
  {
    const auto & face = _grid.face(iface);
    _pmap.insert({face.center(), iface});
  }
  auto it_dest = face_indices.begin();
  for (auto it = _pmap.begin(); it != _pmap.end(); ++it, ++it_dest)
    *it_dest = it->second;
  _pmap.clear();
}

}  // end namespace discretization
