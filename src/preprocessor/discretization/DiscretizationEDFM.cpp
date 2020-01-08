#include "DiscretizationEDFM.hpp"
#include "DiscretizationDFM.hpp"

namespace discretization
{

DiscretizationEDFM::
DiscretizationEDFM(const DoFNumbering & split_dof_numbering,
                   const DoFNumbering & combined_dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data)
    : m_split_dofs(split_dof_numbering),
      DiscretizationBase(combined_dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationEDFM::build()
{
  DiscretizationDFM discr_dfm(m_split_dofs, m_data, m_split_cv, m_split_con);
  discr_dfm.build();
  build_control_volume_data_();
  // build_connection_data_();
}

void DiscretizationEDFM::build_control_volume_data_()
{
  for (auto & cv : m_cv_data)
  {
    cv.volume = 0;
    cv.center = {0.0, 0.0, 0.0};
  }

  // first compute parent volumes since some props are weighted by them
  for (const auto & cv : m_split_cv)
  {
    size_t parent_dof;
    if (cv.type == ControlVolumeType::cell)
      parent_dof = m_dofs.cell_dof(cv.master);
    else // if (cv.type == ControlVolumeType::face)
      parent_dof = m_dofs.face_dof( cv.master );
    m_cv_data[parent_dof].volume += cv.volume;
  }

  for (const auto & cv : m_split_cv)
  {
    size_t parent_dof;
    if (cv.type == ControlVolumeType::cell)
    {
      parent_dof = m_dofs.cell_dof(cv.master);
    }
    else // if (cv.type == ControlVolumeType::face)
      parent_dof = m_dofs.face_dof( cv.master );

    auto & parent_cv = m_cv_data[parent_dof];
    parent_cv.type = cv.type;
    const double volume_fraction = cv.volume / parent_cv.volume;
    parent_cv.center += cv.center * volume_fraction;
    parent_cv.porosity += cv.porosity * volume_fraction;
    parent_cv.permeability = cv.permeability;  // assume they are the same
    parent_cv.custom = cv.custom;  // assume they are the same
  }
}

void DiscretizationEDFM::build_connection_data_()
{
  // if (cv_data.size() == 0) return;
  // const size_t max_edfm_index = m_min_edfm_index + cv_data.size() - 1;
  // for (const auto & con: m_con)
  // {
  //   const auto edfm_elements = find_edfm_elements_(con);
  //   if (con.type == discretization::ConnectionType::matrix_fracture)
  //   {
  //     if (!edfm_elements.empty())
  //       build_matrix_fracture_(con);
  //   }
  //   else if (con.type == discretization::ConnectionType::fracture_fracture)
  //   {
  //     if (edfm_elements.size() == con.elements.size())
  //       build_edfm_edfm_(con);
  //     else if (edfm_elements.empty())
  //       build_edfm_dfm_(con);
  //   }
  // }
  // convert_flow_map_to_vector_();
}

void DiscretizationEDFM::build_matrix_fracture_(const ConnectionData &con)
{
  // assert(con.elements.size() == 2);
  // // first element is fracture
  // // std::cout << con.elements[0] << std::endl;
  // assert(m_cv[con.elements[0]].type == discretization::ControlVolumeType::face);

  // const auto &cell = m_grid.cell(m_cv[con.elements[1]].master);
  // // edfm intersections might split cell in more than 2 sub-cells
  // const auto &parent_cell = cell.ultimate_parent();
  // // connect to cell parent since we split edfm cells
  // size_t con_index;
  // if (m_con_map.contains(con.elements[0], con.elements[1]))
  //   con_index = m_con_map.index(con.elements[0], con.elements[1]);
  // else
  //   con_index = m_con_map.insert(con.elements[0], con.elements[1]);
  // auto &new_con = m_con_map.get_data(con_index);

  // // transmissibility is weighted average of two halves m-F transmissibilities
  // if (new_con.coefficients.empty()) new_con.coefficients = {0.0, 0.0};
  // new_con.coefficients[0] += con.coefficients[0] * cell.volume() / parent_cell.volume();
  // new_con.coefficients[1] -= - new_con.coefficients[0];
}

void DiscretizationEDFM::build_edfm_edfm_(const ConnectionData & con)
{
  assert(con.elements.size() == 2);
  std::cout << "edfm-edfm " << con.elements[0] << " " << con.elements[1] << std::endl;

  size_t con_index;
  if (m_con_map.contains(con.elements[0], con.elements[1]))
    con_index = m_con_map.index(con.elements[0], con.elements[1]);
  else
    con_index = m_con_map.insert(con.elements[0], con.elements[1]);
  auto &new_con = m_con_map.get_data(con_index);
  new_con = con;
}

void DiscretizationEDFM::build_edfm_dfm_(const ConnectionData & con)
{
  // assert (m_min_edfm_index <= con.elements[1] && m_min_edfm_index + m_n_edfm_faces > con.elements[1]);
  // assert(con.elements.size() == 2);
  // std::cout << "edfm-dfm " << con.elements[0] << " " << con.elements[1] << std::endl;

  // // first element is dfm, second is dfm
  // size_t con_index;
  // if (m_con_map.contains(con.elements[0], con.elements[1]))
  //   con_index = m_con_map.index(con.elements[0], con.elements[1]);
  // else
  //   con_index = m_con_map.insert(con.elements[0], con.elements[1]);
  // auto &new_con = m_con_map.get_data(con_index);
  // new_con = con;
}

void DiscretizationEDFM::convert_flow_map_to_vector_()
{
  // mcon_data.reserve( m_con_map.size() );
  // for (auto it = m_con_map.begin(); it != m_con_map.end(); ++it)
  // {
  //   con_data.emplace_back();
  //   auto & con = con_data.back();
  //   con = *it;
  // }
}

}  // end namespace discretization
