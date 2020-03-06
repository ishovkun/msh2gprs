#ifdef WITH_EIGEN
#include "PolyhedralElementDirect.hpp"
#include "EdgeComparison.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "EdgeComparison.hpp"
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
  compute_vertex_mapping_();

  std::cout << "done building discr" << std::endl;
  std::cout << "_element_grid.n_active_cells " << _element_grid.n_active_cells() << std::endl;

  build_face_boundary_conditions_();
}

void PolyhedralElementDirect::compute_vertex_mapping_()
{
  std::vector<size_t> node_tags;
  std::vector<double> node_coord, parametric_coord;
  gmsh::model::mesh::getNodes(node_tags, node_coord, parametric_coord, /*dim = */ 3,
                              /* tag */ -1, /*includeBoundary =*/ true,
                              /* return_parametric =  */ false);
  auto it = std::max_element(node_tags.begin(), node_tags.end());
  assert( it != node_tags.end() );
  const size_t max_node_tag = *it + 1;
  _vertex_mapping.resize(max_node_tag, std::numeric_limits<size_t>::max());
  size_t iv = 0;
  for (const size_t node : node_tags)
    _vertex_mapping[node]  = iv++;
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
    // initialize system matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> face_system_matrix =
        Eigen::SparseMatrix<double,Eigen::RowMajor>(vertex_dofs.size(), vertex_dofs.size());
    // fill system matrix
    build_face_system_matrix_(iface, face_system_matrix, vertex_dofs);
    const std::vector<size_t> parent_face_vertices = _parent_cell.vertices();
    for (size_t ipv=0; ipv<parent_face_vertices.size(); ++ipv)
    {
      Eigen::SparseMatrix<double,Eigen::RowMajor> face_system_matrix_with_bc = face_system_matrix;
      Eigen::VectorXd rhs = Eigen::VectorXd::Zero( vertex_dofs.size() );
      impose_bc_on_face_system_( ipv, vertex_dofs, face_system_matrix_with_bc, rhs );
      face_system_matrix_with_bc.makeCompressed();

      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(face_system_matrix_with_bc);
      solver.factorize(face_system_matrix_with_bc);
      solver.analyzePattern(face_system_matrix_with_bc);
      solver.factorize(face_system_matrix_with_bc);
      solver.solve(rhs);
    }


    exit(0);
  }
}

void PolyhedralElementDirect::impose_bc_on_face_system_(const size_t parent_vertex,
                                                        const std::vector<size_t> & vertex_dofs,
                                                        Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
                                                        Eigen::VectorXd & rhs)
{
  for (size_t iv=0; iv<_support_edge_vertices[parent_vertex].size(); ++iv)
  {
    const size_t vertex = _support_edge_vertices[parent_vertex][iv];
    const size_t dof = vertex_dofs[vertex];
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, dof); it; ++it)
      it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
    rhs[dof] = _support_edge_values[parent_vertex][iv];
  }
}

void PolyhedralElementDirect::build_face_system_matrix_(const size_t iface,
                                                        Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                                        const std::vector<size_t> & vertex_dofs)
{
  const int element_type = api::get_gmsh_element_id(angem::VTK_ID::TriangleID);  // triangle
  // NOTE: gmsh labels faces starting from 1, ergo iface + 1
  const int gmsh_plane_marker = iface + 1;
  FeValues fe_values( element_type, gmsh_plane_marker );

  // get imformation on the elements we're gonna be looping over
  std::vector<size_t> element_node_tags;
  std::vector<size_t> tags;
  gmsh::model::mesh::getElementsByType( element_type, tags, element_node_tags, gmsh_plane_marker );

  const size_t nv = api::get_n_vertices(element_type);
  Eigen::MatrixXd cell_matrix(fe_values.n_vertices(), fe_values.n_vertices());
  std::vector<size_t> face_dofs(nv);

  for (size_t ielement=0; ielement<tags.size(); ++ielement)
  {
    for (size_t iv=0; iv<nv; ++iv)
    {
      const size_t vertex_tag = element_node_tags[nv * ielement + iv];
      const size_t vertex = _vertex_mapping[vertex_tag];
      assert (vertex < _element_grid.n_vertices() );
      face_dofs[iv] = vertex_dofs[ vertex ];
      assert( face_dofs[iv] !=  std::numeric_limits<size_t>::max());
    }

    cell_matrix.setZero();
    fe_values.update(ielement);

    // assemble local element matrix
    for (size_t q = 0; q < fe_values.n_q_points(); ++q)
    {
      for (size_t i = 0; i < nv; ++i)
        for (size_t j = 0; j < nv; ++j)
          cell_matrix(i, j) += -(fe_values.grad(i, q) * // grad phi_i(x_q)
                                 fe_values.grad(j, q) * // grad phi_j(x_q)
                                 fe_values.JxW(q));     // dV
    }

    /* distribute local to global */
    for (size_t i = 0; i < nv; ++i)
      for (size_t j = 0; j < nv; ++j)
      {
        const size_t idof = face_dofs[i];
        const size_t jdof = face_dofs[j];
        // std::cout << "ij = " << idof << " " << jdof << std::endl;
        face_system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
      }
  }
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

void PolyhedralElementDirect::build_edge_boundary_conditions_()
{
  std::vector<std::vector<size_t>> vertex_parent_faces  = map_vertices_to_parent_faces_();
  std::vector<std::list<size_t>> parent_vertex_faces = map_parent_vertices_to_parent_faces_();
  const auto pair_markers_to_edge = edgecmp::EdgeComparison::get_edges( parent_vertex_faces );
  const auto parent_nodes = _parent_cell.polyhedron()->get_points();
  _support_edge_vertices.resize( parent_nodes.size() );
  _support_edge_values.resize( parent_nodes.size() );

  for (size_t v=0; v<_element_grid.n_vertices(); ++v)
  {
    if ( vertex_parent_faces[v].size() > 1 )
    {
      const auto & markers = vertex_parent_faces[v];
      for (size_t mi=0; mi<markers.size(); ++mi)
        for (size_t mj=mi+1; mj<markers.size(); ++mj)
        {
          const size_t m1 = markers[mi];
          const size_t m2 = markers[mj];
          assert( pair_markers_to_edge.contains(m1, m2) );
          const auto & parent_edge = pair_markers_to_edge.get_data(m1, m2);
          const size_t vp1 = parent_edge.either();
          const size_t vp2 = parent_edge.other(vp1);
          const Point p1 = parent_nodes[vp1];
          const Point p2 = parent_nodes[vp2];
          const double dp = p1.distance(p2);
          const double d1 = _element_grid.vertex(v).distance(p1);
          const double d2 = _element_grid.vertex(v).distance(p2);
          _support_edge_vertices[vp1].push_back(v);
          _support_edge_values[vp1].push_back( (dp - d1 ) / dp );
          _support_edge_vertices[vp2].push_back(v);
          _support_edge_values[vp2].push_back( (dp - d2) / dp );

          // also set zero on edges that are not either parent
          for (size_t vp=0; vp<parent_nodes.size(); ++vp)
            if (vp != vp1 && vp != vp2)
            {
              _support_edge_vertices[vp].push_back(v);
              _support_edge_values[vp].push_back(0.0);
            }
        }
    }
  //   // else if (vertex_markers[v].size() == 1)
  //   // {
  //   //   const size_t iface = vertex_markers[v].front();
  //   //   for (size_t pv = 0; pv < _cell.n_vertices(); pv++)
  //   //   {
  //   //     if ( std::find( parent_vertex_markers[pv].begin(), parent_vertex_markers[pv].end(),
  //   //                     iface) == parent_vertex_markers[pv].end() )
  //   //     {
  //   //       _support_edge_vertices[pv].push_back(v);
  //   //       _support_edge_values[pv].push_back(0.0);
  //   //     }
  //   //   }
  //   // }
  }
}

std::vector<std::vector<size_t>> PolyhedralElementDirect::map_vertices_to_parent_faces_()
{
  /* Build structure that store face markers for each vertex in grid */
  std::vector<std::vector<size_t>> vertex_markers( _element_grid.n_vertices() );
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 1 && face->marker() > 0)
    {
      const size_t ipf = face->marker();  // i parent face
      for (size_t v : face->vertices())
        if (std::find( vertex_markers[v].begin(), vertex_markers[v].end(), ipf ) ==
            vertex_markers[v].end())
          vertex_markers[v].push_back(ipf);
    }
  return vertex_markers;
}

std::vector<std::list<size_t>> PolyhedralElementDirect::map_parent_vertices_to_parent_faces_()
{
  // map parent vertex to parent faces
  std::vector<std::list<size_t>> parent_vertex_markers( _parent_cell.n_vertices() );
  const std::vector<size_t> parent_vertices = _parent_cell.sorted_vertices();
  const std::vector<const mesh::Face*> parent_faces = _parent_cell.faces();
  for (size_t ipf=0; ipf<parent_faces.size(); ++ipf)
  {
    const auto parent_face = parent_faces[ipf];

    for (size_t iv_parent=0; iv_parent<parent_vertices.size(); ++iv_parent)
    {
      const size_t parent_vertex = parent_vertices[iv_parent];
      if (parent_face->has_vertex(parent_vertex))
        parent_vertex_markers[iv_parent].push_back( ipf + 1);
    }
  }
  return parent_vertex_markers;
}



}  // end namespace discretization

#endif  // with_eigen
