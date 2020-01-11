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
  if (m_cv_data.size() != m_dofs.n_dofs())
    m_cv_data.resize( m_dofs.n_dofs() );
}

void DiscretizationBase::build_cell_data_()
{
  for (auto cell = m_grid.begin_active_cells(); cell != m_grid.end_active_cells(); ++cell)
  {
    const std::size_t i = cell->index();
    auto & data = m_cv_data[ m_dofs.cell_dof(i) ];
    data.type = ControlVolumeType::cell;
    data.master = i;
    data.porosity = m_data.get_porosity(i);
    data.permeability = m_data.get_permeability(i);
    data.center = cell->center();
    data.volume = cell->volume() * data.porosity;

    data.custom.resize(m_data.output_flow_properties.size());
    for (size_t j = 0; j < m_data.output_flow_properties.size(); ++j)
      data.custom[j] = m_data.cell_properties[m_data.output_flow_properties[j]][i];
  }
}

}
