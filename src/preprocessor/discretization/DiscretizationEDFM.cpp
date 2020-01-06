#include "DiscretizationEDFM.hpp"

namespace discretization
{

DiscretizationEDFM::
DiscretizationEDFM(const std::vector<DiscreteFractureConfig> & dfm_fractures,
                   const std::vector<discretization::ControlVolumeData> & mixed_cv_data,
                   const std::vector<discretization::ConnectionData> & mixed_connection_data,
                   const size_t n_dfm_faces, // const size_t n_cells,
                   gprs_data::SimData & data)
    : DiscretizationBase(dfm_fractures, data),
      m_cv(mixed_cv_data),
      m_con(mixed_connection_data),
      m_min_edfm_index(n_dfm_faces),
      m_n_edfm_faces(calculate_edfm_faces_())
{}

void DiscretizationEDFM::build()
{
  extract_control_volume_data_();
  build_connection_data_();
}

size_t DiscretizationEDFM::calculate_edfm_faces_() const
{
  size_t n_edfm_faces = 0;
  for (std::size_t i=m_min_edfm_index; i<m_cv.size(); ++i)
    if ( m_cv[i].type == discretization::ControlVolumeType::face)
      n_edfm_faces++;
  assert ( m_min_edfm_index + n_edfm_faces > 0 );
  return n_edfm_faces;
}

void DiscretizationEDFM::extract_control_volume_data_()
{
  for (std::size_t i=m_min_edfm_index; i<m_min_edfm_index + m_n_edfm_faces; ++i)
    if ( m_cv[i].type == discretization::ControlVolumeType::face)
      cv_data.push_back( m_cv[i] );
}

void DiscretizationEDFM::build_connection_data_()
{
  if (cv_data.size() == 0) return;
  const size_t max_edfm_index = m_min_edfm_index + cv_data.size() - 1;
  for (const auto & con: m_con)
    if (con.type == discretization::ConnectionType::matrix_fracture)
    {
  // std::cout << "here " << con.elements[0] << " " << con.elements[1] << std::endl;
      // if (m_min_edfm_index <= con.elements[0] && max_edfm_index >= con.elements[0])
      if (!find_edfm_elements_(con).empty())
        build_matrix_fracture_(con, m_min_edfm_index, max_edfm_index);
    }
    else if (con.type == discretization::ConnectionType::fracture_fracture)
      build_fracture_fracture_(con);

  assert ( false && "write extraction code" );
}

void DiscretizationEDFM::build_matrix_fracture_(const ConnectionData &con,
                                                const size_t min_edfm_index,
                                                const size_t max_edfm_index)
{
  assert(con.elements.size() == 2);
  // first element is fracture
  std::cout << con.elements[0] << std::endl;
  assert(m_cv[con.elements[0]].type == discretization::ControlVolumeType::face);

  const auto &cell = m_grid.cell(m_cv[con.elements[1]].master);
  // edfm intersections might split cell in more than 2 sub-cells
  const auto &parent_cell = m_grid.cell(cell.ultimate_parent());
  // connect to cell parent since we split edfm cells
  size_t con_index;
  if (m_con_map.contains(con.elements[0], con.elements[1]))
    con_index = m_con_map.index(con.elements[0], con.elements[1]);
  else
    con_index = m_con_map.insert(con.elements[0], con.elements[1]);
  auto &new_con = m_con_map.get_data(con_index);
  // transmissibility is weighted average of two halves m-F transmissibilities
  new_con.coefficients.resize(2);
  new_con.coefficients[0] =
      con.coefficients[0] * cell.volume() / parent_cell.volume();
  new_con.coefficients[1] = new_con.coefficients[0];
}

std::vector<size_t> DiscretizationEDFM::find_edfm_elements_(const ConnectionData & con)
{
  std::vector<size_t> result;
  for (size_t i = 0; i < con.elements.size(); i++ )
    if (m_min_edfm_index <= con.elements[i] && m_min_edfm_index + m_n_edfm_faces > con.elements[i])
      result.push_back(i);
  return result;
  
}

void DiscretizationEDFM::build_fracture_fracture_(const ConnectionData & con)
{
  // find_edfm_elements_(con);
}

}  // end namespace discretization
