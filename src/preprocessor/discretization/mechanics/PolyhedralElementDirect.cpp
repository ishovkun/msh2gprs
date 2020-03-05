#ifdef WITH_EIGEN
#include "PolyhedralElementDirect.hpp"
#include "EdgeComparison.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "VTKWriter.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

namespace discretization {

using api = gprs_data::GmshInterface;
using FeValues = gprs_data::FeValues;
using Point = angem::Point<3,double>;

PolyhedralElementDirect::PolyhedralElementDirect(const mesh::Cell & cell)
    : _parent_cell(cell)
{
  build_();
}

void PolyhedralElementDirect::build_()
{
  api::build_triangulation(_parent_cell, _element_grid);
  std::cout << "done building discr" << std::endl;
  std::cout << "_element_grid.n_active_cells " << _element_grid.n_active_cells() << std::endl;

  build_face_boundary_conditions_();

}

void PolyhedralElementDirect::build_face_boundary_conditions_()
{
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> face_domains = create_face_domains_();
  const std::vector<const mesh::Face*> parent_faces = _parent_cell.faces();
  const std::vector<size_t> parent_vertices = _parent_cell.vertices();
  for (size_t iface=0; iface<face_domains.size(); ++iface)
  {
    // identify vertices the will constitute the linear system and create dof mapping
    const std::vector<size_t> vertex_dofs = create_face_vertex_dofs_(face_domains[iface]);
    build_face_system_matrix_(iface, vertex_dofs);
    // for (size_t pv = 0; pv < parent_vertices.size(); pv++)
    // {
     
    // }
  }
}

void PolyhedralElementDirect::build_face_system_matrix_(const size_t iface, const std::vector<size_t> & vertex_dofs)
{
  // NOTE: 2 is the triangle gmsh id
  // NOTE: gmsh labels faces starting from 1, ergo iface + 1
  FeValues fe_values( 2, iface + 1 );

}


std::vector<std::vector<size_t>> PolyhedralElementDirect::create_face_domains_()
{
  std::vector<std::vector<size_t>> parent_face_children(_parent_cell.faces().size());
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const size_t parent_face_index = face->marker() - 1;
      parent_face_children[ parent_face_index ].push_back( face->index() );
    }

  return parent_face_children;
}

std::vector<size_t> PolyhedralElementDirect::create_face_vertex_dofs_(const std::vector<size_t> & face_indices)
{
  const size_t unmarked = std::numeric_limits<size_t>::max();
  std::vector<size_t> dofs(_element_grid.n_vertices(), unmarked);
  size_t dof = 0;
  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = _element_grid.face(iface);
    for (const size_t v : face.vertices())
      if ( dofs[v] == unmarked)
        dofs[v] = dof++;
  }

  return dofs;
}

}  // end namespace discretization

#endif  // with_eigen
