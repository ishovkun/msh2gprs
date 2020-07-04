#pragma once

#include "mesh/Mesh.hpp"
#include <map>

namespace discretization {

struct point_comparator
{
  explicit point_comparator(const angem::Plane<double> & pln)
      : _pln(pln) {};

  bool operator() (const angem::Point<3,double> & p1, const angem::Point<3,double> &p2) const
  {
    // const auto pp1 = _pln.project_point(p1);
    // const auto pp2 = _pln.project_point(p2);
    const auto &pp1 = p1;
    const auto &pp2 = p2;
    for (int i = 2; i >=0; i--)
    {
      const double diff = pp1[i] - pp2[i];
      if (std::fabs(diff) < tol)
        continue;
      else if (diff < 0)
        return true;
      else // if (diff > 0)
        return false;
    }
    return false;
  }

  const double tol = 1e-9;
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
