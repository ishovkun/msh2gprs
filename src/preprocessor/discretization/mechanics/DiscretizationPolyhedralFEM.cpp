#include "DiscretizationPolyhedralFEM.hpp"
#ifdef WITH_EIGEN
#include "Elements/PolyhedralElementDirect.hpp"
#include "Elements/PolyhedralElementMSRSB.hpp"
#endif

namespace discretization {

DiscretizationPolyhedralFEM::DiscretizationPolyhedralFEM(mesh::Mesh & grid, const FiniteElementConfig & config,
                                                         const std::vector<int> & fracture_face_markers,
                                                         const std::vector<size_t> & neumann_face_indices)
    : DiscretizationFEMBase::DiscretizationFEMBase(grid, config, fracture_face_markers, neumann_face_indices)
{}

void DiscretizationPolyhedralFEM::build_(mesh::Cell & cell)
{
#ifdef WITH_EIGEN
  if (_config.solver == SolverType::direct || _config.solver == SolverType::cg)
    _element = std::make_shared<PolyhedralElementDirect>(cell, _grid, _config);
  else if (_config.solver == SolverType::msrsb)
    _element = std::make_shared<PolyhedralElementMSRSB>(cell, _grid, _config);
#elseif
  throw std::runtime_error("Eigen is not available");
#endif
}

}  // end namespace discretization
