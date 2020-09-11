#ifdef WITH_EIGEN
#include "PolyhedralElementDirect.hpp"
#include "EdgeComparison.hpp"  // provides Edge and EdgeComparison
#include "mesh/Subdivision.hpp"  // provides mesh::Subdivision
#include "FeValues.hpp"          // provides discretization::FeValues
#include "EdgeComparison.hpp"
#include "VTKWriter.hpp"         // debugging, provides io::VTKWriter
#include "FEMFaceDoFManager.hpp" // provides discretization::FEMFaceDoFManager
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

namespace discretization {

using api = gprs_data::GmshInterface;
using Point = angem::Point<3,double>;
const size_t UNMARKED = std::numeric_limits<size_t>::max();


PolyhedralElementDirect::PolyhedralElementDirect(const mesh::Cell & cell,
                                                 const mesh::Mesh & grid,
                                                 const FiniteElementConfig & config)
    : PolyhedralElementBase(cell, grid, config, true)
{
  // build grid for the element
  build_triangulation_();
  // solve problems on faces
  try {
    build_face_boundary_conditions_();
    // construct the laplace system matrix for the cell volume laplace equation
    build_cell_system_matrix_();
    // impose BC's and solve laplace system to get shape functions
    compute_shape_functions_();
  }
  catch (const std::exception& e)
  {
    std::cout << "Error: " << e.what() << std::endl;
    save_face_domains_(std::to_string(cell.index()) + ".vtk");
    throw e;
  }
}

void PolyhedralElementDirect::build_face_boundary_conditions_()
{
  build_edge_boundary_conditions_();
  // identify child faces that belong to each face parent
  _face_domains = create_face_domains_();
  const std::vector<const mesh::Face*> parent_faces = _parent_cell.faces();
  const std::vector<size_t> parent_vertices = _parent_cell.vertices();
  FEMFaceDoFManager dof_manager;
  for (size_t iface=0; iface<_face_domains.size(); ++iface)
  {
    assert(_face_domains[iface].size() > 2);

    // identify vertices the will constitute the linear system and create dof mapping
    const DoFNumbering vertex_numbering = dof_manager.build(_subgrid, _face_domains[iface]);
    // initialize system matrix
    auto face_system_matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>(vertex_numbering.n_dofs(),
                                                                          vertex_numbering.n_dofs());
    // fill system matrix
    build_face_system_matrix_(iface, face_system_matrix, _face_domains[iface], vertex_numbering);
    const std::vector<size_t> parent_face_vertices = parent_faces[iface]->vertices();
    face_system_matrix.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(face_system_matrix);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero( vertex_numbering.n_dofs() );
    for (size_t ipv=0; ipv<parent_face_vertices.size(); ++ipv)
    {
      rhs.setZero();
      const size_t pv = std::distance(parent_vertices.begin(),
                                      std::find( parent_vertices.begin(), parent_vertices.end(),
                                                 parent_face_vertices[ipv] ));
      if (ipv == 0)
        impose_bc_on_face_system_( pv, vertex_numbering, face_system_matrix, rhs,
                                   /*impose_on_matrix = */ true);
      else
        impose_bc_on_face_system_( pv, vertex_numbering, face_system_matrix, rhs,
                                   /*impose_on_matrix = */ false);

      if (ipv == 0)  // factorize only once since matrix doesn't change, only rhs
      {
        solver.analyzePattern(face_system_matrix);
        solver.factorize(face_system_matrix);
      }
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
        for (size_t vertex = 0; vertex < _subgrid.n_vertices(); ++vertex)
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
                                                        Eigen::VectorXd & rhs,
                                                        const bool impose_on_matrix)
{
  // for (size_t iv=0; iv<_support_edge_vertices[parent_vertex].size(); ++iv)
  // {
  //   const size_t vertex = _support_edge_vertices[parent_vertex][iv];
  //   if (vertex_dofs.has_vertex(vertex))
  //   {
  //     const size_t dof = vertex_dofs.vertex_dof(vertex);
  //     for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, dof); it; ++it)
  //       it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
  //     rhs[dof] = _support_edge_values[parent_vertex][iv];
  //   }

  // }
  for (size_t vertex = 0; vertex < _subgrid.n_vertices(); ++vertex)
    if (vertex_dofs.has_vertex(vertex) && _support_edge_values[parent_vertex][vertex] >= 0)
    {
      const size_t dof = vertex_dofs.vertex_dof(vertex);
      if (impose_on_matrix)
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, dof); it; ++it)
          it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
      rhs[dof] = _support_edge_values[parent_vertex][vertex];
    }
}

void PolyhedralElementDirect::build_cell_system_matrix_()
{
  // build element jacobian for the homogeneous laplace equation
  _system_matrix = Eigen::SparseMatrix<double,Eigen::ColMajor>(_subgrid.n_vertices(),
                                                               _subgrid.n_vertices());
  // since we only build tetrahedral element mesh
  const size_t nv = 4;
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  Eigen::MatrixXd cell_matrix(nv, nv);
  for (auto cell = _subgrid.begin_active_cells(); cell != _subgrid.end_active_cells(); ++cell)
  {
    cell_matrix.setZero();
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

void PolyhedralElementDirect::build_face_system_matrix_(const size_t parent_face,
                                                        Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                                        const std::vector<size_t> & face_indices,
                                                        const DoFNumbering & vertex_dofs)
{
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  const size_t nv = 3;
  Eigen::MatrixXd cell_matrix(nv, nv);
  std::vector<size_t> face_dofs(nv);

  // all faces within the parent face must have the same orientation
  const auto basis = _subgrid.face(face_indices.front()).polygon().plane().get_basis();
  fe_values.set_basis(basis);

  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = _subgrid.face(iface);
    fe_values.update(face);
    cell_matrix.setZero();

    // get vertex dofs
    const auto & verts = face.vertices();
    for (size_t iv=0; iv<nv; ++iv)
      face_dofs[iv] = vertex_dofs.vertex_dof( verts[iv] );

    // assemble local element matrix
    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
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


void PolyhedralElementDirect::build_edge_boundary_conditions_()
{
  std::vector<std::vector<size_t>> vertex_parent_faces  = map_vertices_to_parent_faces_();
  std::vector<std::list<size_t>> parent_vertex_faces = map_parent_vertices_to_parent_faces_();
  const auto pair_markers_to_edge = edgecmp::EdgeComparison::get_edges( parent_vertex_faces,
                                                                        _parent_cell.polyhedron()->get_edges());
  const auto parent_nodes = _parent_cell.polyhedron()->get_points();
  _support_edge_vertices.resize( parent_nodes.size() );
  _support_edge_values.resize( parent_nodes.size(), std::vector<double>(_subgrid.n_vertices(), -1.) );

  for (size_t v=0; v<_subgrid.n_vertices(); ++v)
  {
    if ( vertex_parent_faces[v].size() > 1 )
    {
      const auto & markers = vertex_parent_faces[v];
      for (size_t mi=0; mi<markers.size(); ++mi)
        for (size_t mj=mi+1; mj<markers.size(); ++mj)
        {
          const size_t m1 = markers[mi];
          const size_t m2 = markers[mj];
          if ( !pair_markers_to_edge.contains(m1, m2) )  // happend in some patrially split cells
          {
            continue;
          }
          const auto & edges = pair_markers_to_edge.get_data(m1, m2);
          for (const auto & parent_edge : edges)
          {
            const size_t vp1 = parent_edge.either();
            const size_t vp2 = parent_edge.other(vp1);
            const Point p1 = parent_nodes[vp1];
            const Point p2 = parent_nodes[vp2];
            const double dp = p1.distance(p2);
            const double d1 = _subgrid.vertex(v).distance(p1);
            const double d2 = _subgrid.vertex(v).distance(p2);
            // _support_edge_vertices[vp1].push_back(v);
            // _support_edge_values[vp1].push_back( (dp - d1 ) / dp );
            // _support_edge_vertices[vp2].push_back(v);
            // _support_edge_values[vp2].push_back( (dp - d2) / dp );
            const double val1 = (dp - d1) / dp;
            const double val2 = (dp - d2) / dp;
            // if < 0 it means it's a vertex  on an edge with internal parent vertex
            // (hangin node)
            // then dirichlet value sohuld be zero
            if (val1 < 0 || val2 < 0)
            {
              _support_edge_values[vp1][v] = 0;
              _support_edge_values[vp2][v] = 0;
            }
            else
            {
              _support_edge_values[vp1][v] = val1;
              _support_edge_values[vp2][v] = val2;
            }
            // if ( (dp - d1) / dp < 0 || (dp - d2) / dp < 0)
            // {
            //   std::cout << "suck " << v << " " << vp1 << " " << vp2 << std::endl;
            //   throw "fuck";
            // }
            // assert( (dp - d1) / dp >= 0.0 );
            // assert( (dp - d2) / dp >= 0.0 );

            // if ( vp1 == 2 || vp2 == 2 )
            // if (v == 2 && (vp1 == 2 || vp2 == 2))
            // {
            //   std::cout << "Set phi["<<vp1<<"] = " << (dp - d1 ) / dp
            //             << " in " << v
            //             << " on " << vp1 << "-" << vp2
            //             << std::endl;
            //   std::cout << "Set phi["<<vp2<<"] = " <<  (dp - d2) / dp
            //             << " in " << v
            //             << " on " << vp1 << "-" << vp2
            //             << std::endl;
            // }
            // also set zero on edges that do not have either parent
            for (size_t vp=0; vp<parent_nodes.size(); ++vp)
              if (vp != vp1 && vp != vp2)
              {
                if (std::fabs(_support_edge_values[vp][v] + 1) < 1e-8)
                  _support_edge_values[vp][v] = 0.0;
                // _support_edge_vertices[vp].push_back(v);
                // _support_edge_values[vp].push_back(0.0);
                // if (v == 2 && vp == 2)
                // std::cout << "Set phi["<<vp<<"] = " << 0
                //           << " in " << v
                //           << " on " << vp1 << "-" << vp2
                //           << std::endl;
              }
          }
        }
    }
  }
}

std::vector<std::vector<size_t>> PolyhedralElementDirect::map_vertices_to_parent_faces_()
{
  /* Build structure that store face markers for each vertex in grid */
  std::vector<std::vector<size_t>> vertex_markers( _subgrid.n_vertices() );
  for (auto face = _subgrid.begin_active_faces(); face != _subgrid.end_active_faces(); ++face)
    if (face->neighbors().size() == 1 && face->marker() > 0)
    {
      const size_t ipf = face->marker();  // index parent face
      for (size_t v : face->vertices())
        if (std::find( vertex_markers[v].begin(), vertex_markers[v].end(), ipf ) ==
            vertex_markers[v].end())
          vertex_markers[v].push_back(ipf);
    }
  return vertex_markers;
}

void PolyhedralElementDirect::debug_save_boundary_face_solution(const std::string fname) const
{
  std::cout << "saving " << fname << std::endl;

  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_subgrid, out);
  const size_t nv = _subgrid.n_vertices();
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
  for (size_t i=0; i<_subgrid.n_vertices(); ++i)
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
    _basis_functions.push_back(Eigen::VectorXd::Zero(_subgrid.n_vertices()));
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero( _subgrid.n_vertices() );

  // can compress it since imposing bc's operates on row iterators
  _system_matrix.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(_system_matrix);


  for (size_t pv=0; pv<n_parent_vertices; ++pv)
  {
    rhs.setZero();
    //  impose boundary conditions
    //  NOTE: no need to copy the matrix since BC's locations don't change
    if (pv == 0)
      impose_boundary_conditions_(_system_matrix, rhs, pv);
    else
      impose_boundary_conditions_(rhs, pv);

    // TODO: I cannot currently use conjugate gradient since the matrix is non-symmetric
    // due to dirichlet conditions. I need to symmetrize the matrix by modifying the
    // columns that correspond to dirichlet vertices
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Lower|Eigen::Upper> solver;
    // solver.analyzePattern(mat);
    // solver.factorize(mat);
    // solver.setMaxIterations(200);
    // solver.setTolerance(1e-10);
    // solver.compute(mat);

    if ( pv == 0 )  // need to factorize only once since only rhs changes
    {
      solver.analyzePattern(_system_matrix);
      solver.factorize(_system_matrix);
    }

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

void PolyhedralElementDirect::impose_boundary_conditions_(Eigen::VectorXd & rhs, const size_t ipv)
{
  for (size_t iv=0; iv<_support_boundary_vertices[ipv].size(); ++iv)
  {
    const size_t i = _support_boundary_vertices[ipv][iv];
    rhs[i] = _support_boundary_values[ipv][iv];
  }
}

void PolyhedralElementDirect::save_face_domains_(std::string fname)
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_subgrid, out);
  const size_t nv = _subgrid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  std::vector<double> output(nv, 0);
  for (auto face = _subgrid.begin_faces(); face != _subgrid.end_faces(); ++face)
    if (face->marker() != 0)
    {
      for (auto v : face->vertices())
        output[v] = face->marker();
    }
  IO::VTKWriter::add_data(output, "marker", out);

}


}  // end namespace discretization

#endif  // with_eigen
