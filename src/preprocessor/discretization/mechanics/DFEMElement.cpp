#include "DFEMElement.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "VTKWriter.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace discretization {

using api = gprs_data::GmshInterface;
using FeValues = gprs_data::FeValues;

DFEMElement::DFEMElement(const mesh::Cell & cell)
    : _cell(cell)
{
  build_();
}

void DFEMElement::build_()
{
  build_triangulation_();
  build_jacobian_();
  initial_guess_();
  debug_save_shape_functions_();
  run_msrsb_();
}

void DFEMElement::build_triangulation_()
{
  api::build_triangulation(_cell);
  // ask gmsh to provide gaussian points, shape functions, and jacobians
  api::get_elements(_element_types, _element_tags, _element_nodes, /* dim = */ 3);
  numberNodesEndElements_(_element_types, _element_tags, _element_nodes);

  std::vector<double> crds, prcrds;
  gmsh::model::mesh::getNodes(_node_tags, crds, prcrds, /*dim*/ 3,
                             /* tag */ -1, /*includeBoundary =*/ true, false);
  // convert node coordinates into angem points
  _node_coord.resize(crds.size() / 3);
  for (std::size_t i=0; i<crds.size()/3; ++i)
  {
    angem::Point<3,double> node = { crds[3*i], crds[3*i+1], crds[3*i+2] };
    _node_coord[_node_numbering[_node_tags[i]]] = std::move(node);
  }

  // build support region boundaries for msrsb
  build_support_boundaries_();
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
    gmsh::model::mesh::getNodes(bnd_nodes, crds, prcrds, /*dim*/ 2,
                                /* tag */ parent_face+1,
                                /*includeBoundary =*/ true, false);
    for (std::size_t ivertex=0; ivertex<parent_vertices.size(); ++ivertex)
      if ( faces[parent_face]->has_vertex(parent_vertices[ivertex]) )
        for (const size_t vertex : bnd_nodes)
          _support_boundaries[ivertex].insert( _node_numbering[vertex] );
  }
}

void DFEMElement::run_msrsb_()
{
  Eigen::DiagonalPreconditioner<double> jacobi_preconditioner(_system_matrix);
  // jacobi_preconditioner.setMaxIterations(1);
  std::vector<Eigen::VectorXd> solutions;
  const size_t parent_n_vert = _cell.vertices().size();
  for (std::size_t i=0; i<parent_n_vert; ++i)
    solutions.push_back(Eigen::VectorXd::Zero(_node_numbering.size()));

  for (size_t parent_node = 0; parent_node < parent_n_vert; parent_node++)
  {
    jacobi_preconditioner._solve_impl(_basis_functions[parent_node], solutions[parent_node]);
    // _basis_functions[parent_node] = solutions[i];
  }
   for (size_t vertex=0; vertex < _node_numbering.size(); vertex++)
   {
     enforce_zero_on_boundary_(vertex, solutions);
     enforce_partition_of_unity_(vertex, solutions);
   }

  debug_save_shape_functions_("shape_functions1.vtk");
 
}

void DFEMElement::enforce_partition_of_unity_(const size_t fine_vertex,
                                              std::vector<Eigen::VectorXd> & solutions)
{
  // compute the sum of new values
  double sum_new_values = 0;
  const size_t parent_n_vert = _cell.vertices().size();
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
    sum_new_values += solutions[parent_vertex][fine_vertex];

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

// double DFEMElement::jacobi_iteration_(const size_t fine_vertex)
// {
 
// }

void DFEMElement::build_jacobian_()
{
  // build element jacobian for the homogeneous laplace equation
  _system_matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>(_node_numbering.size(),
                                                              _node_numbering.size());

  for (std::size_t itype = 0; itype < _element_types.size(); ++itype)
  {
    const size_t type = _element_types[itype];
    FeValues fe_values(type, _element_tags[itype].size());

    for (std::size_t itag=0; itag<_element_tags[itype].size(); ++itag)
    {
      const size_t tag = _element_tags[itype][itag];
      const size_t cell_number = _cell_numbering[tag];
      fe_values.update(cell_number);
      Eigen::MatrixXd cell_matrix = Eigen::MatrixXd::Zero(fe_values.n_vertices(),
                                                          fe_values.n_vertices());

      for (size_t q = 0; q < fe_values.n_q_points(); ++q)
      {
        for (size_t i = 0; i < fe_values.n_vertices(); ++i)
          for (size_t j = 0; j < fe_values.n_vertices(); ++j)
            cell_matrix(i, j) += (fe_values.grad(i, q) * // grad phi_i(x_q)
                                  fe_values.grad(j, q) * // grad phi_j(x_q)
                                  fe_values.JxW(q));           // dx
      }

      for (size_t i = 0; i < fe_values.n_vertices(); ++i)
        for (size_t j = 0; j < fe_values.n_vertices(); ++j)
        {
          const size_t inode = _element_nodes[itype][fe_values.n_vertices()*itag + i];
          const size_t jnode = _element_nodes[itype][fe_values.n_vertices()*itag + j];
          const size_t idof = _node_numbering[inode];
          const size_t jdof = _node_numbering[jnode];
          _system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
        }
    }
  }
  _system_matrix.makeCompressed();
}

void DFEMElement::
numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags)
{
  _node_numbering.clear();
  _cell_numbering.clear();
  size_t dof_cell = 0, dof_node = 0;
  for (std::size_t itype = 0; itype < element_types.size(); ++itype)
  {
    const size_t nv = api::get_n_vertices(element_types[itype]);

    for (std::size_t itag=0; itag<element_tags[itype].size(); ++itag)
    {
      _cell_numbering[itag] = dof_cell++;
      for (std::size_t inode=0; inode<nv; ++inode)
      {
        const size_t node = node_tags[itype][nv * itag + inode];
        if (_node_numbering.find(node) == _node_numbering.end())
          _node_numbering[node] = dof_node++;
      }
    }
  }
}

void DFEMElement::initial_guess_()
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (std::size_t i=0; i<parent_n_vert; ++i)
    _basis_functions.push_back(Eigen::VectorXd::Zero(_node_numbering.size()));

  // get nodes that form the vertices of the original elements
  // std::vector<size_t> node_tags;
  // std::vector<double> coord, parametric_coord;
  // gmsh::model::mesh::getNodes(node_tags, coord, parametric_coord,
  //                             0, 0, /* incl boundary */ true,
  //                             /* return parametric = */ false);
  // for (auto & it : node_tags)
  //     std::cout << it << std::endl;
  // std::cout << "coord" << std::endl;
  // for (auto v : coord)
  //   std::cout << v << std::endl;
  // std::cout << "cell v" << std::endl;
  // std::cout << _cell.polyhedron()->get_points()[0] << std::endl;

  // all nodes
  // std::cout << std::endl;
  // std::cout << "all nodes" << std::endl;

  const auto parent_nodes = _cell.polyhedron()->get_points();
  for (std::size_t i=0; i<_node_tags.size(); ++i)
  {
    double max_dist = std::numeric_limits<double>::max();
    size_t closest = 0;
    const angem::Point<3,double> & node = _node_coord[i];
    for (std::size_t j=0; j<parent_nodes.size(); ++j)
    {
      const auto & pn = parent_nodes[j];
      const double dist = pn.distance(node);
      if (dist < max_dist)
      {
        max_dist = dist;
        closest = j;
      }
    }
    _basis_functions[closest][_node_numbering[_node_tags[i]]] = 1.0;
  }
}

void DFEMElement::debug_save_shape_functions_(const std::string fname)
{
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
    std::vector<double> output(_node_coord.size());
    for (std::size_t i=0; i<_node_coord.size(); ++i)
      output[i] = _basis_functions[j][i];
    IO::VTKWriter::add_data(output, "basis-"+std::to_string(j), out);
  }
  out.close();
}

// bool DFEMElement::in_support_boundary_(const size_t fine_vertex, const size_t parent_vertex) const
// {
//   if ( _support_boundaries[parent_face].find(fine_vertex) == _support_boundaries[parent_face].end() )
//     return false;
//   else return true;
// }

// bool DFEMElement::in_global_support_boundary_(const size_t fine_vertex) const
// {
//   for (size_t parent_vertex = 0; parent_vertex < _cell.vertices().size(); parent_vertex++)
//     if (in_support_boundary_(fine_vertex, parent_vertex))
//       return true;
//   return false;
// }

}  // end namespace discretization
