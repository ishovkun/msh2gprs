#pragma once

#include "mesh/Cell.hpp"                    // provices mesh::Cell
#include "FiniteElementData.hpp"            // provides FiniteElementData
#include "FeValues.hpp"


namespace discretization {

/**
 * This class implements a standard finite element discretization for a
 * single cell of the mesh.
 */
class StandardFiniteElement {
 public:
  /**
   * Constructor.
   * Build FEM discretization of the cell.
   * Input:
   * @param[in] cell : grid cell to be discretized
   */
  StandardFiniteElement(const mesh::Cell & cell);
  // get FE data for volume integration
  const FiniteElementData & get_cell_data() const {return _cell_data;}
  // get FE data for surface integration
  const std::vector<FiniteElementData> & get_face_data() const {return _face_data;}

 protected:
  template <angem::VTK_ID vtk_id>
  void build_(FeValues<vtk_id> & fe_values, FiniteElementData & entity_data);

  const mesh::Cell & _cell;
  FiniteElementData _cell_data;                // FEM values and gradients in cell integration points
  std::vector<FiniteElementData> _face_data;   // FEM values and gradients in face integration points
};

template <angem::VTK_ID vtk_id>
void StandardFiniteElement::build_(FeValues<vtk_id> & fe_values, FiniteElementData & entity_data)
{
  const size_t nv = ELEMENT_DIM<vtk_id>;
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
    entity_data.points.push_back(data);

    // center
    FEPointData c;
    c.values.resize( nv, 0.0 );
    c.grads.resize( nv, {0,0,0} );
    for (size_t vertex=0; vertex<nv; ++vertex)
    {
      c.values[vertex] = fe_values.value_center(vertex);
      c.grads[vertex] = fe_values.grad_center(vertex);
    }
    c.weight = fe_values.detJ_center();
    entity_data.center = std::move(c);
  }
 
}


}  // end namespace discretization
