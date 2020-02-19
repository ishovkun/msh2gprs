#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "VTKWriter.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace discretization {

using api = gprs_data::GmshInterface;
using FeValues = gprs_data::FeValues;

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
}

void DFEMElement::build_triangulation_()
{
  api::build_triangulation(_cell);
  // ask gmsh to provide gaussian points, shape functions, and jacobians
  api::get_elements(_element_types, _element_tags, _element_nodes, /* dim = */ 3);

  std::vector<double> crds, prcrds;
  gmsh::model::mesh::getNodes(_node_tags, crds, prcrds, /*dim = */ 3,
                              /* tag */ -1, /*includeBoundary =*/ true,
                              /* return_parametric =  */ false);

  numberNodesEndElements_(_element_types, _element_tags, _element_nodes);
  // convert node coordinates into angem points
  _node_coord.resize(_node_tags.size());
  for (std::size_t i=0; i<_node_tags.size(); ++i)
  {
    angem::Point<3,double> node = { crds[3*i], crds[3*i+1], crds[3*i+2] };
    _node_coord[i] = std::move(node);
  }

  // build support region boundaries for msrsb
  build_support_boundaries_();
  save_support_boundaries_();
}

void DFEMElement::build_support_boundaries_()
{
  const std::vector<const mesh::Face*> faces = _cell.faces();
  const size_t n_parent_faces = faces.size();
  const std::vector<size_t> parent_vertices = _cell.vertices();
  // get nodes for surface
  std::vector<double> crds, prcrds;
  std::vector<size_t> bnd_nodes;
  _support_boundaries.resize(parent_vertices.size());
  for (std::size_t parent_face=0; parent_face<n_parent_faces; ++parent_face)
  {
    bnd_nodes.clear();
    gmsh::model::mesh::getNodes(bnd_nodes, crds, prcrds, /*dim*/ 2,
                                /* tag */ parent_face + 1,
                                /*includeBoundary =*/ true, false);
    for (size_t ivertex=0; ivertex<parent_vertices.size(); ++ivertex)
      if ( !faces[parent_face]->has_vertex(parent_vertices[ivertex]) )
        for (const size_t vertex_tag : bnd_nodes)
          _support_boundaries[ivertex].insert( _node_numbering[vertex_tag] );
  }
}

void DFEMElement::save_support_boundaries_()
{
  const std::string fname = "support_bnd.vtk";
  std::cout << "saving " << fname << std::endl;
  std::vector<std::vector<size_t>> cells;
  for (std::size_t itype = 0; itype < _element_types.size(); ++itype)
  {
    const size_t nv = api::get_n_vertices(_element_types[itype]);
    for (size_t itag = 0; itag < _element_tags[itype].size(); ++itag)
    {
      std::vector<size_t> cell(nv);
      for (std::size_t v=0; v<nv; ++v)
        cell[v] = _node_numbering[_element_nodes[itype][nv * itag + v]];
      cells.push_back( std::move(cell) );
    }
  }

  std::vector<int> cell_types(cells.size(), angem::VTK_ID::TetrahedronID);
  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_node_coord, cells, cell_types, out);
  IO::VTKWriter::enter_section_point_data(_node_coord.size(), out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(_node_coord.size(), 0);
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
    solutions.push_back(Eigen::VectorXd::Zero(_node_numbering.size()));

  double dphi = 10  * _msrsb_tol;
  size_t iter = 0;
  const size_t max_iter = 2000;
  do
  {
    dphi = jacobi_iteration_(solutions, preconditioner);
    if (!(iter % 10))
      std::cout << "iter = " << iter << " dphi = " << dphi << std::endl;
    if (iter++ >= max_iter)
      throw std::runtime_error("msrsb did not converge");
  } while (dphi > _msrsb_tol);

  debug_save_shape_functions_("shape_functions-final.vtk");
}

double DFEMElement::jacobi_iteration_(std::vector<Eigen::VectorXd> & solutions,
                                      const JacobiPreconditioner & prec)
{
  for (size_t parent_node = 0; parent_node < _cell.vertices().size(); parent_node++)
  {
    prec.solve(_system_matrix, _basis_functions[parent_node], solutions[parent_node]);
  }
  for (size_t vertex=0; vertex < _node_numbering.size(); vertex++)
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
  _system_matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>(_node_numbering.size(),
                                                               _node_numbering.size());

  for (std::size_t itype = 0; itype < _element_types.size(); ++itype)
  {
    const size_t type = _element_types[itype];
    // std::cout << "element_type = " << type << std::endl;
    FeValues fe_values(type, _element_tags[itype].size());
    Eigen::MatrixXd cell_matrix(fe_values.n_vertices(), fe_values.n_vertices());

    for (size_t etag = 0; etag < _element_tags[itype].size(); ++etag)
    {
      const size_t tag = _element_tags[itype][etag];
      // std::cout << "element_tag = " << tag << std::endl;
      const size_t cell_number = _cell_numbering[tag];

      fe_values.update(cell_number);
      cell_matrix.setZero();

      for (size_t q = 0; q < fe_values.n_q_points(); ++q)
      {
        for (size_t i = 0; i < fe_values.n_vertices(); ++i)
          for (size_t j = 0; j < fe_values.n_vertices(); ++j)
            cell_matrix(i, j) += -(fe_values.grad(i, q) * // grad phi_i(x_q)
                                   fe_values.grad(j, q) * // grad phi_j(x_q)
                                   fe_values.JxW(q));     // dV
      }

      /* distribute local to global */
      for (size_t i = 0; i < fe_values.n_vertices(); ++i)
        for (size_t j = 0; j < fe_values.n_vertices(); ++j)
        {
          const size_t itag = _element_nodes[itype][fe_values.n_vertices()*etag + i];
          const size_t jtag = _element_nodes[itype][fe_values.n_vertices()*etag + j];
          const size_t idof = _node_numbering[itag];
          const size_t jdof = _node_numbering[jtag];
          _system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
        }
    }
  }
  _system_matrix.makeCompressed();
  // std::cout << _system_matrix << std::endl;
}

void DFEMElement::
numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags)
{
  _node_numbering.clear();
  _cell_numbering.clear();
  size_t dof_cell = 0, dof_node = 0;
  // number elements
  for (std::size_t itype = 0; itype < element_types.size(); ++itype)
    for (std::size_t itag=0; itag<element_tags[itype].size(); ++itag)
      _cell_numbering[itag] = dof_cell++;

  // number nodes
  for (const size_t tag : _node_tags)
    _node_numbering[tag] = dof_node++;
}

void DFEMElement::initial_guess_()
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (std::size_t i=0; i<parent_n_vert; ++i)
    _basis_functions.push_back(Eigen::VectorXd::Zero(_node_numbering.size()));

  const auto parent_nodes = _cell.polyhedron()->get_points();
  for (size_t i=0; i < _node_tags.size(); ++i)
  {
    double min_dist = std::numeric_limits<double>::max();
    size_t closest = 0;
    const angem::Point<3,double> & node = _node_coord[i];
    for (size_t j = 0; j < parent_nodes.size(); ++j)
    {
      const auto & pn = parent_nodes[j];
      const double dist = pn.distance(node);
      if (dist < 1e-9)
      {
        std::cout << "point should not be in support region "<< j << " vertex "<< i <<"(" <<_node_tags[i] << ")" << std::endl;
        assert(!in_support_boundary_(i, j));
      }
      if (dist < min_dist)
      {
        min_dist = dist;
        closest = j;
      }
    }
    _basis_functions[closest][i] = 1.0;
  }
}

void DFEMElement::debug_save_shape_functions_(const std::string fname)
{
  std::cout << "saving " << fname << std::endl;
  std::vector<std::vector<size_t>> cells;
  for (std::size_t itype = 0; itype < _element_types.size(); ++itype)
  {
    const size_t nv = api::get_n_vertices(_element_types[itype]);
    for (std::size_t itag = 0; itag < _element_tags[itype].size(); ++itag) {
      std::vector<size_t> cell(nv);
      for (std::size_t v=0; v<nv; ++v)
        cell[v] = _node_numbering[_element_nodes[itype][nv * itag + v]];
      cells.push_back( std::move(cell) );
    }
  }

  std::vector<int> cell_types(cells.size(), angem::VTK_ID::TetrahedronID);

  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_node_coord, cells, cell_types, out);
  IO::VTKWriter::enter_section_point_data(_node_coord.size(), out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(_node_coord.size(), 0.0);
    for (std::size_t i=0; i<_node_coord.size(); ++i)
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

}  // end namespace discretization

#endif  // with_eigen
