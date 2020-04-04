#ifdef WITH_EIGEN
#include "PolyhedralElementDirect.hpp"
#include "EdgeComparison.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "FeValues.hpp"
#include "EdgeComparison.hpp"
#include "VTKWriter.hpp"
#include "FEMFaceDoFManager.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

namespace discretization {

using api = gprs_data::GmshInterface;
using Point = angem::Point<3,double>;
const size_t UNMARKED = std::numeric_limits<size_t>::max();


PolyhedralElementDirect::PolyhedralElementDirect(const mesh::Cell & cell)
    : _parent_cell(cell)
{
  build_();
}

void PolyhedralElementDirect::build_()
{
  // triangulate the polyhedral element
  api::build_triangulation(_parent_cell, _element_grid);
  // get the relation between gmsh vertex ids and grid vertices
  // compute_vertex_mapping_();
  // solve problems on faces
  build_face_boundary_conditions_();
  // construct the laplace system matrix for the cell volume laplace equation
  build_cell_system_matrix_();
  // impose BC's and solve laplace system to get shape functions
  compute_shape_functions_();
  // TODO: debugging, delete this
  // debug_save_shape_functions_("shape_functions-final.vtk");
  // identify the locations of the gauss points for the polyhedral element
  find_integration_points_();
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  compute_cell_fe_quantities_();
  // compute shape function values, gradients, and weights in the
  // integration points in faces
  _face_data.resize( _face_gauss_points.size() );
  const size_t nfaces = _parent_cell.faces().size();
  for (size_t iface=0; iface<nfaces; ++iface)
    compute_face_fe_quantities_(iface);
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
  _vertex_mapping.resize(max_node_tag, UNMARKED);
  size_t iv = 0;
  for (const size_t node : node_tags)
    _vertex_mapping[node]  = iv++;
}

void PolyhedralElementDirect::build_face_boundary_conditions_()
{
  build_edge_boundary_conditions_();
  // identify child faces that belong to each face parent
  _face_domains = create_face_domains_();
  const std::vector<const mesh::Face*> parent_faces = _parent_cell.faces();
  const std::vector<size_t> parent_vertices = _parent_cell.sorted_vertices();
  FEMFaceDoFManager dof_manager;
  for (size_t iface=0; iface<_face_domains.size(); ++iface)
  {
    // identify vertices the will constitute the linear system and create dof mapping
    const DoFNumbering vertex_numbering = dof_manager.build(_element_grid, _face_domains[iface]);
    // initialize system matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> face_system_matrix =
        Eigen::SparseMatrix<double,Eigen::RowMajor>(vertex_numbering.n_dofs(), vertex_numbering.n_dofs());
    // fill system matrix
    build_face_system_matrix_(iface, face_system_matrix, _face_domains[iface], vertex_numbering);
    const std::vector<size_t> parent_face_vertices = parent_faces[iface]->vertices();
    for (size_t ipv=0; ipv<parent_face_vertices.size(); ++ipv)
    {
      // copy matrix, create rhs vector, and impose bc's on them
      Eigen::SparseMatrix<double,Eigen::RowMajor> face_system_matrix_with_bc = face_system_matrix;
      Eigen::VectorXd rhs = Eigen::VectorXd::Zero( vertex_numbering.n_dofs() );
      const size_t pv = std::distance(parent_vertices.begin(),
                                      std::find( parent_vertices.begin(), parent_vertices.end(),
                                                 parent_face_vertices[ipv] ));
      impose_bc_on_face_system_( pv, vertex_numbering, face_system_matrix_with_bc, rhs );
      face_system_matrix_with_bc.makeCompressed();

      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(face_system_matrix_with_bc);
      solver.analyzePattern(face_system_matrix_with_bc);
      solver.factorize(face_system_matrix_with_bc);
      const Eigen::VectorXd solution = solver.solve(rhs);
      // add the solution to the face BC container
      append_face_solution_(pv, solution, vertex_numbering);
    }

    // append zeroes in face vertices to boundary conditions if face doesn't contain a vertex
    for (size_t pv=0; pv<parent_vertices.size(); ++pv)
    {
      const size_t parent_grid_vertex = parent_vertices[pv];
      if (std::find(parent_face_vertices.begin(), parent_face_vertices.end(),
                    parent_grid_vertex) == parent_face_vertices.end())
      {
        for (size_t vertex = 0; vertex < _element_grid.n_vertices(); ++vertex)
          if (vertex_numbering.has_vertex(vertex)) {
            _support_boundary_vertices[pv].push_back(vertex);
            _support_boundary_values[pv].push_back(0.0);
          }
      }
    }

  }
}

void PolyhedralElementDirect::impose_bc_on_face_system_(const size_t parent_vertex,
                                                        const DoFNumbering & vertex_dofs,
                                                        Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
                                                        Eigen::VectorXd & rhs)
{
  for (size_t iv=0; iv<_support_edge_vertices[parent_vertex].size(); ++iv)
  {
    const size_t vertex = _support_edge_vertices[parent_vertex][iv];
    if (vertex_dofs.has_vertex(vertex))
    {
      const size_t dof = vertex_dofs.vertex_dof(vertex);
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, dof); it; ++it)
        it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
      rhs[dof] = _support_edge_values[parent_vertex][iv];
    }

  }
}

void PolyhedralElementDirect::build_cell_system_matrix_()
{
  // build element jacobian for the homogeneous laplace equation
  _system_matrix = Eigen::SparseMatrix<double,Eigen::ColMajor>(_element_grid.n_vertices(),
                                                               _element_grid.n_vertices());

  // since we only build tetrahedral element mesh
  const int element_type = api::get_gmsh_element_id(angem::VTK_ID::TetrahedronID);
  const size_t nv = 4;
  FeValues<angem::VTK_ID::TetrahedronID> fe_values(_element_grid);

  Eigen::MatrixXd cell_matrix(nv, nv);
  for (auto cell = _element_grid.begin_active_cells(); cell != _element_grid.end_active_cells(); ++cell)
  {
    cell_matrix.setZero();
    // fe_values.update(cell->index());
    fe_values.update(*cell);

    const std::vector<size_t> & cell_vertices = cell->vertices();
    const size_t nv = cell_vertices.size();

    // This assembles a local matrix for the laplace equation
    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
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
        const size_t idof = cell_vertices[i];
        const size_t jdof = cell_vertices[j];
        _system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
      }
  }
}

void PolyhedralElementDirect::build_face_system_matrix_(const size_t iface,
                                                        Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                                        const DoFNumbering & vertex_dofs)
{
  // face discretizaiton is triangles
  const int element_type = api::get_gmsh_element_id(angem::VTK_ID::TriangleID);
  // NOTE: gmsh labels faces starting from 1, ergo iface + 1
  const int gmsh_plane_marker = iface + 1;
  gprs_data::FeValues fe_values( element_type, gmsh_plane_marker );

  // get imformation on the elements we're gonna be looping over
  std::vector<size_t> element_node_tags;
  std::vector<size_t> tags;
  gmsh::model::mesh::getElementsByType( element_type, tags, element_node_tags, gmsh_plane_marker );

  const size_t nv = api::get_n_vertices(element_type);
  Eigen::MatrixXd cell_matrix(nv, nv);
  std::vector<size_t> face_dofs(nv);
  for (size_t ielement=0; ielement<tags.size(); ++ielement)
  {
    for (size_t iv=0; iv<nv; ++iv)
    {
      const size_t vertex_tag = element_node_tags[nv * ielement + iv];
      const size_t vertex = _vertex_mapping[vertex_tag];
      assert (vertex < _element_grid.n_vertices() );
      face_dofs[iv] = vertex_dofs.vertex_dof( vertex );
      assert(vertex_dofs.has_vertex(vertex));
    }

    cell_matrix.setZero();
    fe_values.update(ielement);

    // assemble local element matrix
    for (size_t q = 0; q < fe_values.n_q_points(); ++q)
    {
      for (size_t i = 0; i < nv; ++i)
        for (size_t j = 0; j < nv; ++j)
          cell_matrix(i, j) += (fe_values.grad(i, q) * // grad phi_i(x_q)
                                fe_values.grad(j, q) * // grad phi_j(x_q)
                                fe_values.JxW(q));     // dV
    }

    /* distribute local to global */
    for (size_t i = 0; i < nv; ++i)
      for (size_t j = 0; j < nv; ++j)
      {
        const size_t idof = face_dofs[i];
        const size_t jdof = face_dofs[j];
        face_system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
      }
  }
}

void PolyhedralElementDirect::build_face_system_matrix_(const size_t parent_face,
                                                        Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                                        const std::vector<size_t> & face_indices,
                                                        const DoFNumbering & vertex_dofs)
{
  FeValues<angem::VTK_ID::TriangleID> fe_values(_element_grid);
  const size_t nv = 3;
  Eigen::MatrixXd cell_matrix(nv, nv);
  std::vector<size_t> face_dofs(nv);

  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = _element_grid.face(iface);
    fe_values.update(face);
    cell_matrix.setZero();
    // assemble local element matrix
    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
    {
      for (size_t i = 0; i < nv; ++i)
        for (size_t j = 0; j < nv; ++j)
          cell_matrix(i, j) += (fe_values.grad(i, q) * // grad phi_i(x_q)
                                fe_values.grad(j, q) * // grad phi_j(x_q)
                                fe_values.JxW(q));     // dV
      // get vertex dofs
      const auto & verts = face.vertices();
      for (size_t iv=0; iv<nv; ++iv)
        face_dofs[iv] = vertex_dofs.vertex_dof( verts[iv] );

      /* distribute local to global */
      for (size_t i = 0; i < nv; ++i)
        for (size_t j = 0; j < nv; ++j)
        {
          const size_t idof = face_dofs[i];
          const size_t jdof = face_dofs[j];
          face_system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
        }
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

          // also set zero on edges that do not have either parent
          for (size_t vp=0; vp<parent_nodes.size(); ++vp)
            if (vp != vp1 && vp != vp2)
            {
              _support_edge_vertices[vp].push_back(v);
              _support_edge_values[vp].push_back(0.0);
            }
        }
    }
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

void PolyhedralElementDirect::debug_save_boundary_face_solution(const std::string fname) const
{
  std::cout << "saving " << fname << std::endl;

  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _parent_cell.vertices().size();

  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0);
    for (size_t iv=0; iv<_support_boundary_vertices[j].size(); ++iv)
    {
      const size_t v = _support_boundary_vertices[j][iv];
      output[v] = _support_boundary_values[j][iv];
    }

    IO::VTKWriter::add_data(output, "support-" + std::to_string(j), out);
  }
  out.close();
}

void PolyhedralElementDirect::append_face_solution_(const size_t pv, const Eigen::VectorXd & solution,
                                                    const DoFNumbering & vertex_numbering)
{
  _support_boundary_vertices.resize( _parent_cell.vertices().size() );
  _support_boundary_values.resize( _parent_cell.vertices().size() );

  _support_boundary_vertices[pv].reserve( _support_boundary_vertices[pv].size() + vertex_numbering.n_dofs() );
  _support_boundary_values[pv].reserve( _support_boundary_vertices[pv].size() + vertex_numbering.n_dofs() );
  for (size_t i=0; i<_element_grid.n_vertices(); ++i)
    if ( vertex_numbering.has_vertex(i) )
    {
      _support_boundary_vertices[pv].push_back( i );
      _support_boundary_values[pv].push_back( solution[ vertex_numbering.vertex_dof(i) ] );
    }
}

void PolyhedralElementDirect::compute_shape_functions_()
{
  const size_t n_parent_vertices = _parent_cell.vertices().size();
  // allocate shape functions
  for (std::size_t i=0; i<n_parent_vertices; ++i)
    _basis_functions.push_back(Eigen::VectorXd::Zero(_element_grid.n_vertices()));

  for (size_t pv=0; pv<n_parent_vertices; ++pv)
  {
    // assemble problem with BC's for the parent vertex
    // copy system matrix and create a rhs vector
    Eigen::SparseMatrix<double,Eigen::RowMajor> mat = _system_matrix;
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero( _element_grid.n_vertices() );
    // impose boundary conditions
    impose_boundary_conditions_(mat, rhs, pv);
    mat.makeCompressed();

    // TODO: I cannot currently use conjugate gradient since the matrix is non-symmetric
    // due to dirichlet conditions. I need to symmetrize the matrix by modifying the
    // columns that correspond to dirichlet vertices
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Lower|Eigen::Upper> solver;
    // solver.analyzePattern(mat);
    // solver.factorize(mat);
    // solver.setMaxIterations(200);
    // solver.setTolerance(1e-10);
    // solver.compute(mat);

    // Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::RowMajor>> solver;
    // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(mat);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(mat);
    // solver.factorize(mat);
    solver.analyzePattern(mat);
    solver.factorize(mat);

    _basis_functions[pv] = solver.solve(rhs);
    // if ( solver.info() ==  Eigen::ComputationInfo::NumericalIssue)
    //   std::cout << "numerical issue" << std::endl;
    // else if( solver.info() ==  Eigen::ComputationInfo::NoConvergence)
    //   std::cout << "no convergence" << std::endl;
    if (solver.info() != Eigen::ComputationInfo::Success)
    {
      throw std::runtime_error( "Conjugate gratient did not converge" );
    }
  }
}

void PolyhedralElementDirect::impose_boundary_conditions_(Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
                                                          Eigen::VectorXd & rhs, const size_t ipv)
{
  for (size_t iv=0; iv<_support_boundary_vertices[ipv].size(); ++iv)
  {
    const size_t i = _support_boundary_vertices[ipv][iv];
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, i); it; ++it)
      it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
    rhs[i] = _support_boundary_values[ipv][iv];
  }
}

void PolyhedralElementDirect::debug_save_shape_functions_(const std::string fname) const
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _parent_cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0.0);
    for (size_t i = 0; i < nv; ++i)
      output[i] = _basis_functions[j][i];
    IO::VTKWriter::add_data(output, "basis-"+std::to_string(j), out);
  }
  out.close();
}

void PolyhedralElementDirect::compute_cell_fe_quantities_()
{
  // allocate vector of point data
  const size_t n_parents = _basis_functions.size();
  // const int element_type = api::get_gmsh_element_id(angem::VTK_ID::TetrahedronID);
  // gprs_data::FeValues fe_values( element_type );
  FeValues<angem::VTK_ID::TetrahedronID> fe_values(_element_grid);
  for (auto cell = _element_grid.begin_active_cells(); cell != _element_grid.end_active_cells(); ++cell)
    for (size_t q=0; q<_cell_gauss_points.size(); ++q)  //
      if ( cell->polyhedron()->point_inside(_cell_gauss_points[q]) )
      {
        std::vector<angem::Point<3,double>> points = {_cell_gauss_points[q]};
        const std::vector<size_t> & cell_verts = cell->vertices();
        // fe_values.update(cell->index(), points);
        fe_values.update(*cell, _cell_gauss_points[q]);
        const size_t nv = cell->vertices().size();
        const size_t q_loc = 0;  // only one integration point
        FEPointData data;
        data.values.resize( _basis_functions.size() );
        data.grads.resize( _basis_functions.size() );
        for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
        {
          data.values[parent_vertex] = 0;
          for (size_t vertex=0; vertex<nv; ++vertex)
          {
            data.values[parent_vertex] += fe_values.value( vertex, q_loc ) *
                                          _basis_functions[parent_vertex][cell_verts[vertex]];
            data.grads[parent_vertex] += fe_values.grad(vertex, q_loc) *
                                         _basis_functions[parent_vertex][cell_verts[vertex]];
          }
        }

        if ( q < _cell_gauss_points.size() - 1 )
        {
          data.weight = 1.0 / static_cast<double>( _cell_gauss_points.size() - 1 ); // -1 since center is not a gauss point
          _cell_data.points.push_back( std::move(data) );
        }
        else // last point is center
        {
          data.weight = 1.0;
          _cell_data.center = std::move(data);
        }
      }
}

void PolyhedralElementDirect::find_integration_points_()
{
  /* Split a parent cell into pyramids and find their centers */
  const auto polyhedron = _parent_cell.polyhedron();
  std::vector<Point> parent_vertices = polyhedron->get_points();
  parent_vertices.push_back(_parent_cell.center());
  for (const auto & face : polyhedron->get_faces())
    _cell_gauss_points.push_back(create_pyramid_and_compute_center_(face, parent_vertices));
  // add center point to the end of gauss list
  _cell_gauss_points.push_back( _parent_cell.center() );

  /* Loop parent faces, split each face into triangles, and find their centers */
  _face_gauss_points.clear();
  for (const auto & face_poly : polyhedron->get_face_polygons())
  {
    std::vector<angem::Point<3,double>> face_gauss_points =
        split_into_triangles_and_compute_center_(face_poly);
    _face_gauss_points.push_back( face_gauss_points );
    // add face center point
    _face_gauss_points.back().push_back(face_poly.center());
  }
}

Point PolyhedralElementDirect::create_pyramid_and_compute_center_(const std::vector<size_t> & face,
                                                                  const std::vector<Point> & vertices) const
{
  const size_t vertex_center = vertices.size() - 1;
  std::vector<std::vector<size_t>> pyramid_faces;
  for (size_t iv=0; iv<face.size(); ++iv)
  {
    size_t v1, v2;
    v1 = face[iv];
    if ( iv < face.size() - 1 )
      v2 = face[iv+1];
    else
      v2 = face[0];
    pyramid_faces.push_back( {v1, v2, vertex_center} );
  }

  angem::Polyhedron<double> pyramid(vertices, pyramid_faces);
  return pyramid.center();
}

std::vector<angem::Point<3,double>>
PolyhedralElementDirect::split_into_triangles_and_compute_center_(const angem::Polygon<double> & poly)
{
  const angem::Point<3,double> center = poly.center();
  const std::vector<angem::Point<3,double>> & vertices = poly.get_points();
  std::vector<angem::Point<3,double>> result;
  for (size_t i=0; i<vertices.size(); ++i)
  {
    size_t v1, v2;
    if ( i <  vertices.size() - 1)
    {
      v2 = i + 1;
    }
    else
    {
      v2 = 0;
    }

    v1 = i;

    std::vector<angem::Point<3,double>> triangle_vertices =
        {vertices[v1], vertices[v2], center};
    angem::Polygon<double> triangle(triangle_vertices);
    result.push_back( triangle.center() );
  }
  return result;
}

void PolyhedralElementDirect::compute_face_fe_quantities_(const size_t parent_face)
{
  const std::vector<const mesh::Face*> parent_faces = _parent_cell.faces();
  const std::vector<size_t> parent_face_vertices = parent_faces[parent_face]->vertices();
  const std::vector<size_t> parent_vertices = _parent_cell.vertices();
  FeValues<angem::VTK_ID::TriangleID> fe_values(_element_grid);
  const size_t nv = 3; // n vertices in triangle
  const size_t q_loc = 0;  // only one integration point
  for (const size_t iface : _face_domains[parent_face])
  {
    const mesh::Face & face = _element_grid.face(iface);
    const angem::Polygon<double> face_poly = face.polygon();
    const std::vector<size_t> & face_verts = face.vertices();

    const size_t nq = _face_gauss_points[parent_face].size();
    for (size_t q=0; q<nq; ++q)
      if (face_poly.point_inside(_face_gauss_points[parent_face][q]))
      {
        const Point coord = _face_gauss_points[parent_face][q];
        fe_values.update(face, coord);
        FEPointData data;
        data.values.resize( parent_face_vertices.size(), 0.0 );
        data.grads.resize( parent_face_vertices.size(), {0,0,0} );
        for (size_t ipv=0; ipv<parent_face_vertices.size(); ++ipv)
        {
          const size_t pv = std::distance(parent_vertices.begin(),
                                          std::find( parent_vertices.begin(), parent_vertices.end(),
                                                     parent_face_vertices[ipv] ));
          for (std::size_t vertex=0; vertex<nv; ++vertex)
          {
            data.values[ipv] += fe_values.value( vertex, q_loc ) *
                                _basis_functions[pv][face_verts[vertex]];
            data.grads[ipv] += fe_values.grad(vertex, q_loc) *
                               _basis_functions[pv][face_verts[vertex]];
          }
        }
        if ( q < nq - 1 )
        {
          data.weight = 1.0 / static_cast<double>( nq - 1 );
          _face_data[parent_face].points.push_back( std::move(data) );
        }
        else
        {
          data.weight = 1.0;
          _face_data[parent_face].center = std::move(data);
        }
      }

  }


}

}  // end namespace discretization

#endif  // with_eigen
