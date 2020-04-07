#include "StandardFiniteElement.hpp"
#include "angem/VTK_ID.hpp"
#include "FeValues.hpp"

namespace discretization {

using VTK_ID = angem::VTK_ID;

StandardFiniteElement::StandardFiniteElement(const mesh::Cell & cell)
    : _cell(cell)
{
  build_();
}

void StandardFiniteElement::build_()
{
  const VTK_ID id = static_cast<VTK_ID>(_cell.vtk_id());
  switch (id)
  {
    case VTK_ID::TetrahedronID:
      {
        FeValues<VTK_ID::TetrahedronID> fe_values;
        const size_t nv = 4;
        fe_values.update(_cell);
        for (size_t q=0; q<fe_values.n_integration_points(); ++q)
        {
          FEPointData data;
          data.values.resize( nv, 0.0 );
          data.grads.resize( nv, {0,0,0} );
          for (size_t vertex=0; vertex<nv; ++vertex)
          {
            data.values[vertex] = fe_values.value(vertex, q);
            data.grads[vertex] = fe_values.grad(vertex, q);
          }
          data.weight = fe_values.JxW(q);

          // center

          _cell_data.points.push_back(data);

        }

      }
      break;

    case VTK_ID::HexahedronID:
      {
        FeValues<VTK_ID::HexahedronID> fe_values;
        const size_t nv = 6;
        fe_values.update(_cell);
      }
      break;

    default:
      throw std::invalid_argument("FEM for this vtk id is not implemented");
      break;
  }
}

}  // end namespace discretization
