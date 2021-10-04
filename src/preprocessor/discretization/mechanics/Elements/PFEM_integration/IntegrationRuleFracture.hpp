#pragma once
#include "IntegrationRule2d.hpp"

#ifdef WITH_EIGEN

namespace discretization {

class IntegrationRuleFracture
{
 public:
  IntegrationRuleFracture(PolyhedralElementBase const & element,
                          TributaryRegion2dBase const & tributary,
                          size_t                        parent_face,
                          angem::Basis<3, double> const & basis);

  FiniteElementData integrate(std::vector<angem::Point<3,double>> const & vertices,
                              angem::Basis<3, double>             const & basis) const;

  virtual ~IntegrationRuleFracture() = default;

 protected:
  void build_region_(PolyhedralElementBase const & element,
                     std::vector<size_t> const & faces,
                     const angem::Basis<3, double> & basis,
                     size_t offset, size_t region);

  void build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                            FEPointData const & master,
                            FEPointData & target,
                            const angem::Basis<3, double> & basis) const;
  void compute_parent_vertices_(mesh::Cell const & cell, size_t parent_face);
  // returns detJ for face transformation
  double compute_face_scaling_(std::vector<angem::Point<3,double>> const & ref_grad,
                               std::vector<angem::Point<3,double>> const & vertex_coord,
                               angem::Basis<3, double>             const & basis) const;
  double compute_inverse_cell_jacobian_(std::vector<angem::Point<3,double>> const & ref_grad,
                                        angem::Tensor2<3, double>                 & du_dx,
                                        std::vector<angem::Point<3,double>> const & vertex_coord) const;
  void scale_data_(FEPointData & data) const;

  size_t _nregions;
  size_t _npcv;  // number of parent cell vertices
  size_t _npfv;  // number of parent face vertices
  std::vector<size_t> _parent_vertices;  // indices of face vertices in cell::vertices
  std::vector<size_t> _region;
  std::vector<FEPointData> _data;
};

}  // end namespace discretization

#endif
