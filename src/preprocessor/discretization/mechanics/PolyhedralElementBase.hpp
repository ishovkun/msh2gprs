#pragma once
#include "FiniteElementBase.hpp"            // provides FiniteElementBase
#include "config/FiniteElementConfig.hpp"   // provides FiniteElementConfig
#include "mesh/Mesh.hpp"                    // provides mesh::Mesh, mesh::cell
#include "gmsh_interface/GmshInterface.hpp" // provides GmshInterface
#ifdef WITH_EIGEN
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd

namespace discretization {

class IntegrationRuleFacesAverage;
class TributaryRegion2dBase;

/**
 * This is a base class for all polyhedral finite elements.
 * It only performs triangulation
 */
class PolyhedralElementBase : public FiniteElementBase
{
 public:
  const mesh::Mesh & get_grid() const { return _subgrid; }
  // Save a vtk file with shape function values
  void save_shape_functions(const std::string fname) const;
  // Compute cell 3d integration data and return it
  virtual FiniteElementData get_cell_data() override
  {
    if (_cell_data.points.empty())
      build_fe_cell_data_();
    return _cell_data;
  }

  // get FE data for surface integration
  virtual FiniteElementData get_face_data(const size_t iface,
                                          const angem::Basis<3,double> basis) override;
  // get FE data of cell shape functions at face integration points.
  // This is needed for modeling discrete fractures
  virtual FiniteElementData get_fracture_data(const size_t iface,
                                              const angem::Basis<3,double> basis) override;

 protected:
  PolyhedralElementBase(const mesh::Cell & cell,
                        const mesh::Mesh & parent_grid,
                        const FiniteElementConfig & config,
                        const bool sort_faces = false);
  void build_triangulation_();
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> create_face_domains_();
  // map vertices of parent cell to the markers of parent cell face
  std::vector<std::list<size_t>> map_parent_vertices_to_parent_faces_();
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  void build_fe_cell_data_();
  // compute shape function values, gradients, and weights in the
  // integration points in faces
  void build_fe_face_data_();
  // compute shape function values, gradients, and weights in the
  // integration points in fractures
  void build_fe_fracture_data_();
  //
  void build_tributary_2d_(const size_t parent_face);

  const mesh::Cell & _parent_cell;                             // reference to the discretized cell
  const mesh::Mesh & _parent_grid;                             // grid the discrefized cell belongs to
  const FiniteElementConfig & _config;
  const bool _sort_faces;
  mesh::Mesh _subgrid;                                    // triangulation of the discretized cell
  std::vector<Eigen::VectorXd> _basis_functions;               // numerical shape function values
  std::vector<std::vector<size_t>> _face_domains;              // child face indices for each parent face
  std::vector<angem::Point<3,double>> _cell_gauss_points;      // FEM gauss points
  std::vector<std::vector<angem::Point<3,double>>> _face_gauss_points; // FEM face gauss points
  // std::vector<std::vector<angem::Polygon<double>>> _tributary_2d;      // face tributary regions
  std::vector<std::shared_ptr<TributaryRegion2dBase>> _tributary_2d;

  friend class TributaryRegion3dFaces;
  friend class TributaryRegion3dVertices;
  friend class TributaryRegion2dFaces;
  friend class TributaryRegion2dVertices;
  friend class IntegrationRule3dAverage;
  friend class IntegrationRule2dBase;
  friend class IntegrationRule2dAverage;
  friend class IntegrationRule2dPointwise;
  friend class IntegrationRule2dFull;
  friend class IntegrationRuleFractureAverage;
  friend class IntegrationRule3dPointwise;
  friend class IntegrationRuleFractureFull;
  friend class IntegrationRule3dFull;
};

}  // end namespace discretization

#endif
