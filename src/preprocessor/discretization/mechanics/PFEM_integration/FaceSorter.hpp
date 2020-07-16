#pragma once

#include "mesh/Mesh.hpp"
#include <map>

namespace discretization {

struct point_comparator
{
  explicit point_comparator(const angem::Plane<double> & pln,
                            const double tolerance = 1e-9)
      : _pln(pln), _tol(tolerance) {};

  bool operator() (const angem::Point<3,double> & p1, const angem::Point<3,double> & p2) const
  {
    for (int i = 2; i >=0; i--)
    {
      const double diff = p1[i] - p2[i];
      if (std::fabs(diff) < _tol)
        continue;
      else if (diff < 0)
        return true;
      else // if (diff > 0)
        return false;
    }
    return false;
  }

 private:
  const double _tol;
  const angem::Plane<double> _pln;
};

class FaceSorter {
 public:
  FaceSorter(const mesh::Face & parent_face, const mesh::Mesh & subgrid) noexcept;
  virtual ~FaceSorter() = default;
  void sort(std::vector<size_t> & face_indices);

 protected:
  const mesh::Mesh & _grid;
  std::map<angem::Point<3,double>,size_t,point_comparator> _pmap;
};

}  // end namespace discretization
