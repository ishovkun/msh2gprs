#include "DiscretizationStandardFEM.hpp"
#include "Elements/StandardFiniteElement.hpp"

namespace discretization {

DiscretizationStandardFEM::DiscretizationStandardFEM(mesh::Mesh & grid,
                                                     const FiniteElementConfig & config,
                                                     const std::vector<int> & fracture_face_markers,
                                                     const std::vector<int> & neumann_face_markers)
    : DiscretizationFEMBase::DiscretizationFEMBase(grid, config,
                                                   fracture_face_markers,
                                                   neumann_face_markers)
{}

void DiscretizationStandardFEM::build_(mesh::Cell & cell)
{
  _element = std::make_unique<StandardFiniteElement>(cell);
}

}  // end namespace discretization
