#include "FiniteElementBase.hpp"

namespace discretization {

angem::Basis<3,double>
FiniteElementBase::get_face_basis_(mesh::Face const & face,
                                   mesh::Cell const & cell)
{
  angem::Plane plane(face.vertex_coordinates());
  auto fc = face.center();
  auto cc = cell.center();
  if (plane.normal().dot( fc - cc ) < 0)
    plane.get_basis().invert();
  return plane.get_basis();
}

}  // end namespace discretization
