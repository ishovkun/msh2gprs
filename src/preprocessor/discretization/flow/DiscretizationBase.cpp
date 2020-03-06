#include "DiscretizationBase.hpp"

namespace discretization
{

DiscretizationBase::
DiscretizationBase(const DoFNumbering & dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data)
    : m_grid(data.grid),
      m_data(data),
      m_dofs(dof_numbering),
      m_cv_data(cv_data),
      m_con_data(connection_data)
{
  assert ( &m_cv_data );
  if (m_cv_data.size() < m_dofs.n_dofs())
    m_cv_data.resize( m_dofs.n_dofs() );
}

void DiscretizationBase::build_cell_data_(const mesh::Cell& cell)
{
    const std::size_t cell_index = cell.index();
    auto & cv = m_cv_data[ m_dofs.cell_dof(cell_index) ];
    cv.type = ControlVolumeType::cell;
    cv.master = cell_index;
    cv.porosity = m_data.get_porosity(cell_index);
    cv.permeability = m_data.get_permeability(cell_index);
    cv.center = cell.center();
    cv.volume = cell.volume();

    cv.custom.resize(m_data.output_flow_properties.size());
    for (size_t j = 0; j < m_data.output_flow_properties.size(); ++j)
      cv.custom[j] = m_data.cell_properties[m_data.output_flow_properties[j]][cell_index];
}

}
