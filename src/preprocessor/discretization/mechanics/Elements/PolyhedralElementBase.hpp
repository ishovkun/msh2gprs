#pragma once
#include "FiniteElementBase.hpp"            // provides FiniteElementBase
#include "config/FiniteElementConfig.hpp"   // provides FiniteElementConfig
#include "mesh/Mesh.hpp"                    // provides mesh::Mesh, mesh::cell
#include "gmsh_interface/GmshInterface.hpp" // provides GmshInterface
#ifdef WITH_EIGEN
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd

namespace discretization {

class IntegrationRule3d;
class IntegrationRule2d;
class IntegrationRuleFracture;
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
  FiniteElementData get_face_data(size_t iface) override;
  // get FE data of cell shape functions at face integration points.
  // This is needed for modeling discrete fractures
  virtual FiniteElementData get_fracture_data(const size_t iface,
                                              const angem::Basis<3,double> basis) override;
  // enforce cell 3d FE data
  void set_face_cell_data(FiniteElementData const & data) {_cell_data = data;}
  // return the polohedron for the parent cell
  std::unique_ptr<angem::Polyhedron<double>> host_topology() const {return _parent_cell.polyhedron();}
  // returns the reference to the host cell
  mesh::Cell const & host_cell() const { return _parent_cell; }
  // returns raw vectors of basis functions
  std::vector<Eigen::VectorXd> const & get_basis_functions() const noexcept {return _basis_functions;}
  // returns the number of vertices in the host cell
  size_t n_vertices() const {return _parent_cell.n_vertices();}
  // returns the integration scheme
  IntegrationRule3d const & integration_rule3() const {return *_integration_rule3d;}
  // returns reference to the current config
  FiniteElementConfig const & get_config() const {return _config;}
  // get map: parent cell faces to subgrid faces
  std::vector<std::vector<size_t>> const & get_face_domains();
  // get integration rule for cell face iface
  IntegrationRule2d const & integration_rule2(size_t face_local_index);
  // {return *_integration_rules2d[face_local_index];}
  // get integration rule for cell face iface
  IntegrationRuleFracture const & integration_rule_frac(size_t face_local_index);
  // {return *_integration_rules_frac[face_local_index];}

 protected:
  PolyhedralElementBase(const mesh::Cell & cell,
                        const mesh::Mesh & parent_grid,
                        const FiniteElementConfig & config,
                        const bool sort_faces = false);
  // build subgrid of the cell element
  void build_triangulation_();
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> create_face_domains_();
  // map vertices of parent cell to the markers of parent cell face
  std::vector<std::list<size_t>> map_parent_vertices_to_parent_faces_();
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  virtual void build_fe_cell_data_();
  // compute shape function values, gradients, and weights in the
  // integration points in faces
  void build_fe_face_data_();
  // compute shape function values, gradients, and weights in the
  // integration points in fractures
  void build_fe_fracture_data_();
  // build face tributary regions
  std::unique_ptr<TributaryRegion2dBase> build_tributary_2d_(const size_t parent_face);
  // get integration rule used to compute cell fe values

  const mesh::Cell & _parent_cell;                                     // reference to the discretized cell
  const mesh::Mesh & _parent_grid;                                     // grid the discrefized cell belongs to
  const FiniteElementConfig & _config;                                 // configuration
  const bool _sort_faces;                                              // not sure if I need that any more
  mesh::Mesh _subgrid;                                                 // triangulation of the discretized cell
  std::vector<Eigen::VectorXd> _basis_functions;                       // numerical shape function values
  std::vector<std::vector<size_t>> _face_domains;                      // child face indices for each parent face
  std::vector<angem::Point<3,double>> _cell_gauss_points;              // FEM gauss points
  std::vector<std::vector<angem::Point<3,double>>> _face_gauss_points; // FEM face gauss points
  std::shared_ptr<IntegrationRule3d> _integration_rule3d;
  std::vector<std::shared_ptr<IntegrationRule2d>> _integration_rules2d;
  std::vector<std::shared_ptr<IntegrationRuleFracture>> _integration_rules_frac;
};

}  // end namespace discretization

#endif