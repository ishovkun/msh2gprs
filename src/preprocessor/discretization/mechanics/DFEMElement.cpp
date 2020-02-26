#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "VTKWriter.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace discretization {

using api = gprs_data::GmshInterface;
using FeValues = gprs_data::FeValues;
using Point = angem::Point<3,double>;

DFEMElement::DFEMElement(const mesh::Cell & cell,
                         const double msrsb_tol)
    : _cell(cell), _msrsb_tol(msrsb_tol)
{
  build_();
}

void DFEMElement::build_()
{
  build_triangulation_();
  build_jacobian_();
  initial_guess_();
  run_msrsb_();
  debug_save_shape_functions_("shape_functions-final.vtk");
  // save_support_boundaries_();
  find_integration_points_();
  compute_fe_quantities_();
}

void DFEMElement::build_triangulation_()
{
  api::build_triangulation(_cell, _element_grid);
  build_support_boundaries_();
}

void DFEMElement::build_support_boundaries_()
{
  const std::vector<const mesh::Face*> parent_faces = _cell.faces();
  const std::vector<size_t> parent_vertices = _cell.vertices();
  _support_boundaries.resize(parent_vertices.size());
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
  {
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const mesh::Face * parent_face = parent_faces[face->marker() - 1];
      std::vector<size_t> non_adj;  // parent vertices not in parent face
      for (size_t iv_parent=0; iv_parent<parent_vertices.size(); ++iv_parent)
      {
        const size_t parent_vertex = parent_vertices[iv_parent];
        if ( !parent_face->has_vertex(parent_vertex) )
          non_adj.push_back(iv_parent);
      }
      for (const size_t vertex : face->vertices())
        for (const size_t parent_vertex : non_adj)
          _support_boundaries[parent_vertex].insert(vertex);
    }
  }
}

void DFEMElement::save_support_boundaries_()
{
  const std::string fname = "support_bnd.vtk";
  std::cout << "saving " << fname << std::endl;

  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0);
    for (const size_t v : _support_boundaries[j])
      output[v] = 1;
    IO::VTKWriter::add_data(output, "support-"+std::to_string(j), out);
  }
  out.close();
}

void DFEMElement::run_msrsb_()
{
  JacobiPreconditioner preconditioner(_system_matrix);
  std::vector<Eigen::VectorXd> solutions;
  for (size_t i=0; i < _cell.vertices().size(); ++i)
    solutions.push_back(Eigen::VectorXd::Zero( _element_grid.n_vertices() ));

  double dphi = 10  * _msrsb_tol;
  size_t iter = 0;
  const size_t max_iter = 2000;
  do
  {
    dphi = jacobi_iteration_(solutions, preconditioner);

    if (!(iter % 10))  // li'l progress prinout
      std::cout << "iter = " << iter << " dphi = " << dphi << std::endl;

    if (iter++ >= max_iter)
      throw std::runtime_error("msrsb did not converge");

  } while (dphi > _msrsb_tol);
}

double DFEMElement::jacobi_iteration_(std::vector<Eigen::VectorXd> & solutions,
                                      const JacobiPreconditioner & prec)
{
  for (size_t parent_node = 0; parent_node < _cell.vertices().size(); parent_node++)
  {
    prec.solve(_system_matrix, _basis_functions[parent_node], solutions[parent_node]);
  }
  for (size_t vertex=0; vertex < _element_grid.n_vertices(); vertex++)
  {
    enforce_zero_on_boundary_(vertex, solutions);
    enforce_partition_of_unity_(vertex, solutions);
  }

  double error = 0;
  for (size_t parent_node = 0; parent_node < _cell.vertices().size(); parent_node++)
  {
    _basis_functions[parent_node] += solutions[parent_node];
    const double value = solutions[parent_node].norm();
    error = std::max(value, error);
  }
  return error;
}

void DFEMElement::enforce_partition_of_unity_(const size_t fine_vertex,
                                              std::vector<Eigen::VectorXd> & solutions)
{
  const size_t parent_n_vert = _cell.vertices().size();
  double sum_shape_functions = 0;
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
  {
    const double shape = _basis_functions[parent_vertex][fine_vertex];
    const double new_value = shape + solutions[parent_vertex][fine_vertex];
    sum_shape_functions += new_value;
  }
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
  {
    const double shape = _basis_functions[parent_vertex][fine_vertex];
    const double new_value = (shape + solutions[parent_vertex][fine_vertex]) / sum_shape_functions;
    solutions[parent_vertex][fine_vertex] = new_value - shape;
  }
}

void DFEMElement::enforce_zero_on_boundary_(const size_t fine_vertex,
                                            std::vector<Eigen::VectorXd> & solutions)
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
  {
    if (in_support_boundary_(fine_vertex, parent_vertex))
      solutions[parent_vertex][fine_vertex] = 0;
  }
}

void DFEMElement::build_jacobian_()
{
  // build element jacobian for the homogeneous laplace equation
  _system_matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>(_element_grid.n_vertices(),
                                                               _element_grid.n_vertices());
  // since we only build tetrahedral element mesh
  // 4 is tetrahedron id
  // TODO: write a wrapper for it
  FeValues fe_values( 4, _element_grid.n_cells() );

  Eigen::MatrixXd cell_matrix(fe_values.n_vertices(), fe_values.n_vertices());
  for (auto cell = _element_grid.begin_active_cells(); cell != _element_grid.end_active_cells(); ++cell)
  {
    cell_matrix.setZero();
    fe_values.update(cell->index());

    const std::vector<size_t> & cell_vertices = cell->vertices();
    const size_t nv = cell_vertices.size();

    // This assembles a local matrix for the laplace equation
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
        const size_t idof = cell_vertices[i];
        const size_t jdof = cell_vertices[j];
        _system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
      }
  }

  _system_matrix.makeCompressed();
}

void DFEMElement::initial_guess_()
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (std::size_t i=0; i<parent_n_vert; ++i)
    _basis_functions.push_back(Eigen::VectorXd::Zero(_element_grid.n_vertices()));

  const auto parent_nodes = _cell.polyhedron()->get_points();
  for (size_t v=0; v < _element_grid.n_vertices(); ++v)
  {
    double min_dist = std::numeric_limits<double>::max();
    size_t closest = 0;
    for (size_t j = 0; j < parent_nodes.size(); ++j)
    {
      const auto & pn = parent_nodes[j];
      const double dist = pn.distance(_element_grid.vertex(v));
      if (dist < 1e-9)
      {
        std::cout << "point should not be in support region "<< j << " vertex " << v << std::endl;
        assert(!in_support_boundary_(v, j));
      }
      if (dist < min_dist)
      {
        min_dist = dist;
        closest = j;
      }
    }
    _basis_functions[closest][v] = 1.0;
  }
}

void DFEMElement::debug_save_shape_functions_(const std::string fname)
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0.0);
    for (size_t i = 0; i < nv; ++i)
      output[i] = _basis_functions[j][i];
    IO::VTKWriter::add_data(output, "basis-"+std::to_string(j), out);
  }
  out.close();
}

bool DFEMElement::in_support_boundary_(const size_t fine_vertex, const size_t parent_vertex) const
{
  if ( _support_boundaries[parent_vertex].find(fine_vertex) == _support_boundaries[parent_vertex].end() )
    return false;
  else return true;
}

bool DFEMElement::in_global_support_boundary_(const size_t fine_vertex) const
{
  for (size_t parent_vertex = 0; parent_vertex < _cell.vertices().size(); parent_vertex++)
    if (in_support_boundary_(fine_vertex, parent_vertex))
      return true;
  return false;
}


void DFEMElement::find_integration_points_()
{
  /* Split a parent cell into pyramids and find their centers */
  const auto polyhedron = _cell.polyhedron();
  std::vector<Point> parent_vertices = polyhedron->get_points();
  parent_vertices.push_back(_cell.center());
  for (const auto & face : polyhedron->get_faces())
    _integration_points.push_back(create_pyramid_and_compute_center_(face, parent_vertices));
}

Point DFEMElement::create_pyramid_and_compute_center_(const std::vector<size_t> & face,
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

void DFEMElement::compute_fe_quantities_()
{
  // allocate vector of point data
  const size_t n_parents = _basis_functions.size();
  _cell_data.points.resize(_integration_points.size());
  // 4 is a gmsh tetra
  FeValues fe_values( 4, _element_grid.n_cells() );
  for (auto cell = _element_grid.begin_active_cells(); cell != _element_grid.end_active_cells(); ++cell)
    for (size_t q=0; q<_integration_points.size(); ++q)
      if ( cell->polyhedron()->point_inside(_integration_points[q]) )
      {
        std::vector<angem::Point<3,double>> points = {_integration_points[q]};
        const std::vector<size_t> & cell_verts = cell->vertices();
        fe_values.update(cell->index(), points);
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
        data.weight = 1.0 / static_cast<double>( _integration_points.size() );
        _cell_data.points.push_back( std::move(data) );
      }
}

const FiniteElementData & DFEMElement::get_cell_data() const
{
  return _cell_data;
}

}  // end namespace discretization

#endif  // with_eigen
