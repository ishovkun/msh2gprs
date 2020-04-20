#pragma once
#include "FiniteElementBase.hpp"            // provides FiniteElementBase
#include "config/FiniteElementConfig.hpp"   // provides FiniteElementConfig
#include "mesh/Mesh.hpp"                    // provides mesh::Mesh, mesh::cell
#include "gmsh_interface/GmshInterface.hpp" // provides GmshInterface
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd

namespace discretization {

class IntegrationRuleFacesAverage;

/**
 * This is a base class for all polyhedral finite elements.
 * It only performs triangulation
 */
class PolyhedralElementBase : public FiniteElementBase
{
 public:
  const mesh::Mesh & get_grid() const { return _element_grid; }

 protected:
  PolyhedralElementBase(const mesh::Cell & cell, const FiniteElementConfig & config);
  void build_triangulation_();
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> create_face_domains_();

  const mesh::Cell & _parent_cell;                             // reference to the discretized cell
  const FiniteElementConfig & _config;
  mesh::Mesh _element_grid;                                    // triangulation of the discretized cell
  std::vector<Eigen::VectorXd> _basis_functions;               // numerical shape function values
  std::vector<std::vector<size_t>> _face_domains;              // child face indices for each parent face
  std::vector<angem::Point<3,double>> _cell_gauss_points;      // FEM gauss points
  std::vector<std::vector<angem::Point<3,double>>> _face_gauss_points; // FEM face gauss points


  friend class IntegrationRuleFacesAverage;
};

}  // end namespace discretization
