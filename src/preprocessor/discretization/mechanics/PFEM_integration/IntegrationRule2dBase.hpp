#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

class IntegrationRule2dBase {
 public:
  IntegrationRule2dBase(PolyhedralElementBase  & element,
                        const size_t parent_face,
                        const angem::Basis<3, double> & basis)
      : _element(element), _parent_face(parent_face), _basis(basis)
  {
    if (_element._face_domains.empty())
      _element._face_domains = _element.create_face_domains_();
    compute_parent_vertices_();
  }

  virtual ~IntegrationRule2dBase() = default;

  virtual FiniteElementData get() const = 0;

 protected:
  size_t compute_n_parent_vertices_() const noexcept
  {
    return _parent_faces[_parent_face]->vertices().size();
  }

  // map parent face vertices to parent cell vertices
  // to get the index of the basis function
  void compute_parent_vertices_()
  {
    _parent_faces = _element._parent_cell.faces();
    const size_t npv = compute_n_parent_vertices_();
    _parent_vertices.resize(npv);
    const auto parent_cell_vertices = _element._parent_cell.vertices();
    const auto face = _element._parent_cell.faces()[_parent_face];
    for (size_t v=0; v < npv; ++v)
      _parent_vertices[v] =
          std::distance(parent_cell_vertices.begin(),
                        std::find( parent_cell_vertices.begin(), parent_cell_vertices.end(),
                                   face->vertices()[v]));
  }


  PolyhedralElementBase & _element;
  const size_t _parent_face;
  const angem::Basis<3,double> _basis;
  std::vector<const mesh::Face*> _parent_faces;
  std::vector<size_t> _parent_vertices;
};

}  // end namespace discretization


#endif
