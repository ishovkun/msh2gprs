#pragma once

#include "mesh/Cell.hpp"                    // provices mesh::Cell
#include "FiniteElementBase.hpp"            // provides FiniteElementBase
#include "FeValues.hpp"


namespace discretization {

/**
 * This class implements a standard finite element discretization for a
 * single cell of the mesh.
 */
class StandardFiniteElement : public FiniteElementBase {
 public:
  /**
   * Constructor.
   * Build FEM discretization of the cell.
   *
   * The parameter update_face_values should be used if the user needs data for
   * face integration, e.g. Neumann boundary conditions.
   *
   * The parameter udpate_fracture_values dictates whether the user needs
   * to additionally compute cell shape function data at face integration points.
   * These data is used in ADGPRS for DFM fractures, in order to solve the contact problem.
   * Input:
   * @param[in] cell : grid cell to be discretized
   */
  StandardFiniteElement(const mesh::Cell & cell);
  // get FE data for surface integration
  virtual FiniteElementData get_face_data(const size_t iface,
                                          const angem::Basis<3,double> basis) override;
  // get FE data of cell shape functions at face integration points.
  // This is needed for modeling discrete fractures
  virtual FiniteElementData get_fracture_data(const size_t iface,
                                          const angem::Basis<3,double> basis) override;
  virtual ~StandardFiniteElement() = default;

 protected:
  template <angem::VTK_ID vtk_id>
  void build_(FeValues<vtk_id> & fe_values, FiniteElementData & entity_data,
              const bool update_center = true);

  void update_face_values_();
  void update_cell_values_in_faces_();

  const mesh::Cell & _cell;
};

template <angem::VTK_ID vtk_id>
void StandardFiniteElement::build_(FeValues<vtk_id> & fe_values, FiniteElementData & entity_data,
                                   const bool update_center)
{
  const size_t nv = ElementTraits<vtk_id>::n_vertices;
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

    if (update_center)
    {
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
 
}


}  // end namespace discretization
