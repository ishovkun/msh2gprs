#include "DiscretizationBase.hpp"

namespace discretization
{

DiscretizationBase::
DiscretizationBase(const std::vector<DiscreteFractureConfig> & dfm_fractures,
                   gprs_data::SimData & data)
    : m_grid(data.grid),
      m_dfm_config(dfm_fractures),
      m_data(data)
{
  build_dfm_markers_();
}

void DiscretizationBase::build_dfm_markers_()
{
  for (const auto & frac : m_dfm_config)
    m_dfm_markers.insert(frac.label);
  std::cout << "dfm markers: " << std::endl;
  for (int m : m_dfm_markers)
    std::cout << m << " ";
  std::cout << std::endl;
}

bool DiscretizationBase::is_fracture(const int marker) const
{
  const auto it = m_dfm_markers.find(marker);
  if (it != m_dfm_markers.end())
    return true;
  else return false;
}

void DiscretizationBase::build_cell_data_()
{
  cv_data.resize(m_grid.n_cells());
  for (auto cell = m_grid.begin_active_cells(); cell != m_grid.end_active_cells(); ++cell)
  {
    const std::size_t i = cell->index();
    auto & data = cv_data[i];
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


std::vector<ControlVolumeData> & DiscretizationBase::get_cell_data()
{
  return cv_data;
}


std::vector<ConnectionData> & DiscretizationBase::get_face_data()
{
  return con_data;
}


size_t DiscretizationBase::count_dfm_faces_() const
{
  size_t counter = 0;
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (is_fracture( face->marker() ))
      counter++;
  return counter;
}

size_t DiscretizationBase::find_max_cv_index_() const
{
  size_t cv_index = 0;
  // dfm fracs
  for (auto & it_face : m_data.dfm_faces)
    cv_index = std::max( cv_index, it_face.second.cv_index );
  // cells
  for (const size_t cv : m_data.cell_cv_indices)
    cv_index = std::max( cv_index, cv );
  return cv_index;
}

}