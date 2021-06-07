#include "DiscretizationBase.hpp"

namespace discretization
{

DiscretizationBase::
DiscretizationBase(const DoFNumbering & dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data)
    : m_grid(data.grid)
    , m_data(data)
    , m_dofs(dof_numbering)
    , m_cv_data(cv_data)
    , m_con_data(connection_data)
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
    cv.porosity = m_data.cell_properties[m_data.flow.porosity_idx][cell_index];
    auto const &perm_idx = m_data.flow.permeability_idx;
    assert( perm_idx.size() == 3 );
    cv.permeability(0,0) = m_data.cell_properties[ perm_idx[0] ][cell_index];
    cv.permeability(1,1) = m_data.cell_properties[ perm_idx[1] ][cell_index];
    cv.permeability(2,2) = m_data.cell_properties[ perm_idx[2] ][cell_index];
    cv.center = cell.center();
    cv.volume = cell.volume() * m_data.cell_properties[ m_data.flow.vmult_idx ][cell_index];

    cv.custom.resize(m_data.flow.custom_idx.size());
    for (size_t i = 0; i < m_data.flow.custom_idx.size(); ++i)
      cv.custom[i] = m_data.cell_properties[m_data.flow.custom_idx[i]][cell_index];
}

}
