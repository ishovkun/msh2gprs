#pragma once
#include "FiniteElementBase.hpp"            // provides FiniteElementBase
#include "config/FiniteElementConfig.hpp"   // provides FiniteElementConfig
#include "mesh/Mesh.hpp"                    // provides mesh::Mesh, mesh::cell
#include "gmsh_interface/GmshInterface.hpp" // provides GmshInterface
#ifdef WITH_EIGEN
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
  // Save a vtk file with shape function values
  void save_shape_functions(const std::string fname) const;

 protected:
  PolyhedralElementBase(const mesh::Cell & cell,
                        const mesh::Mesh & parent_grid,
                        const FiniteElementConfig & config);
  void build_triangulation_();
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> create_face_domains_();
  // map vertices of parent cell to the markers of parent cell face
  std::vector<std::list<size_t>> map_parent_vertices_to_parent_faces_();


  const mesh::Cell & _parent_cell;                             // reference to the discretized cell
  const mesh::Mesh & _parent_grid;                             // grid the discrefized cell belongs to
  const FiniteElementConfig & _config;
  mesh::Mesh _element_grid;                                    // triangulation of the discretized cell
  std::vector<Eigen::VectorXd> _basis_functions;               // numerical shape function values
  std::vector<std::vector<size_t>> _face_domains;              // child face indices for each parent face
  std::vector<angem::Point<3,double>> _cell_gauss_points;      // FEM gauss points
  std::vector<std::vector<angem::Point<3,double>>> _face_gauss_points; // FEM face gauss points

  friend class TributaryRegion3dFaces;
  friend class TributaryRegion3dVertices;
  friend class TributaryRegion2dFaces;
  friend class TributaryRegion2dVertices;
  friend class IntegrationRule3dAverage;
  friend class IntegrationRule2dAverage;
  friend class IntegrationRule2dPointwise;
  friend class IntegrationRuleFractureAverage;
  friend class IntegrationRule3dPointwise;
  friend class IntegrationRuleFractureMS;
  friend class IntegrationRule3dMS;
};

}  // end namespace discretization

#endif
