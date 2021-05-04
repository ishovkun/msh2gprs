#pragma once
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion3dBase.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

class IntegrationRule3d {
 public:
  IntegrationRule3d(const PolyhedralElementBase & element,
                    const TributaryRegion3dBase  & tributary);
  virtual ~IntegrationRule3d() = default;
  FiniteElementData integrate(std::vector<angem::Point<3,double>> const & vertices) const;

protected:
  void build_region_(PolyhedralElementBase const & element,
                     std::vector<size_t> const & cells,
                     size_t offset,
                     size_t region);
  void build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                            FEPointData const & master,
                            FEPointData & target,
                            angem::Tensor2<3, double> & du_dx) const;
  double compute_detJ_and_invert_cell_jacobian_(const std::vector<angem::Point<3,double>> & ref_grad,
                                                angem::Tensor2<3, double> & du_dx,
                                                std::vector<angem::Point<3,double>> const & vertex_coord) const;
  void scale_data_(FEPointData & data) const;

  size_t _nregions;
  // store full integration data
  std::vector<size_t> _region;
  std::vector<FEPointData> _data;
};

}  // end namespace discretization
