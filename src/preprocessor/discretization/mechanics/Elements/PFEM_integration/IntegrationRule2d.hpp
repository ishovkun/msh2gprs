#pragma once
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion2dBase.hpp"
#include "angem/Tensor2.hpp"

#ifdef WITH_EIGEN

namespace discretization {

class IntegrationRule2d {
 public:
  IntegrationRule2d(const PolyhedralElementBase & element,
                    const TributaryRegion2dBase  & tributary,
                    const size_t parent_face,
                    const angem::Basis<3, double> & basis);

  virtual ~IntegrationRule2d() = default;

  FiniteElementData integrate(std::vector<angem::Point<3,double>> const & vertices,
                              angem::Basis<3, double>             const & basis) const;

 protected:
  void build_region_(PolyhedralElementBase const & element,
                     std::vector<size_t> const & faces,
                     const angem::Basis<3, double> & basis,
                     size_t offset, size_t region);

  void build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                            FEPointData const & master,
                            FEPointData & target,
                            angem::Basis<3, double> const & basis,
                            angem::Tensor2<3, double> & du_dx) const;

  double compute_inverse_jacobian_(std::vector<angem::Point<3,double>> const & ref_grad,
                                   std::vector<angem::Point<3,double>> const & vertex_coord,
                                   angem::Basis<3, double>             const & basis,
                                   angem::Tensor2<3, double> & J_inv) const;
  void scale_data_(FEPointData & data) const;
  void compute_parent_vertices_(mesh::Cell const & cell, size_t parent_face);

  size_t _nregions;
  size_t _npv;  // number of parent vertices
  std::vector<size_t> _parent_vertices;  // indices of face vertices in cell::vertices
  std::vector<size_t> _region;
  std::vector<FEPointData> _data;
};

}  // end namespace discretization

#endif
