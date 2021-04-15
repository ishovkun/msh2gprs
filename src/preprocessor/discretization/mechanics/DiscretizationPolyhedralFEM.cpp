#include "DiscretizationPolyhedralFEM.hpp"
#ifdef WITH_EIGEN
#include "Elements/PolyhedralElementDirect.hpp"
#include "Elements/PolyhedralElementMSRSB.hpp"
#endif

namespace discretization {

DiscretizationPolyhedralFEM::DiscretizationPolyhedralFEM(mesh::Mesh & grid,
                                                         const FiniteElementConfig & config,
                                                         const std::vector<int> & fracture_face_markers,
                                                         const std::vector<int> & neumann_face_markers)
    : DiscretizationFEMBase::DiscretizationFEMBase(grid, config,
                                                   fracture_face_markers,
                                                   neumann_face_markers)
{
}

void DiscretizationPolyhedralFEM::build_(mesh::Cell & cell)
{
  if (_config.solver == SolverType::direct || _config.solver == SolverType::cg)
    _element = std::make_unique<PolyhedralElementDirect>(cell, _grid, _config);
  else if (_config.solver == SolverType::msrsb)
    _element = std::make_unique<PolyhedralElementMSRSB>(cell, _grid, _config);
}

}  // end namespace discretization
