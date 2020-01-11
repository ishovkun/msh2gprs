#include "DiscretizationEDFM.hpp"
#include "DiscretizationDFM.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;

DiscretizationEDFM::
DiscretizationEDFM(const DoFNumbering & split_dof_numbering,
                   const DoFNumbering & combined_dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data,
                   const std::vector<int> & edfm_markers)
    : m_split_dofs(split_dof_numbering),
      m_edfm_markers( edfm_markers.begin(), edfm_markers.end() ),
      DiscretizationBase(combined_dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationEDFM::build()
{
  DiscretizationDFM discr_dfm(m_split_dofs, m_data, m_split_cv, m_split_con);
  discr_dfm.build();
  identify_edfm_faces_();
  build_control_volume_data_();
  build_connection_data_();
}

void DiscretizationEDFM::build_control_volume_data_()
{
  for (auto & cv : m_cv_data)
  {
    cv.volume = 0;
    cv.center = {0.0, 0.0, 0.0};
    cv.aperture = 0;
  }

  // first compute parent volumes since some props are weighted by them
  m_dof_mapping.resize(m_split_cv.size());
  for (size_t i = 0; i < m_split_cv.size(); i++)
    if (m_split_cv[i].volume > 0)  // skip inactive cells
    {
      const auto &cv = m_split_cv[i];
      size_t parent_dof;
      if (cv.type == ControlVolumeType::cell)
        parent_dof = m_dofs.cell_dof(cv.master);
      else // if (cv.type == ControlVolumeType::face)
        parent_dof = m_dofs.face_dof(cv.master);
      m_dof_mapping[i] = parent_dof;
      m_cv_data[parent_dof].volume += cv.volume;
    }

  for (size_t i = 0; i < m_split_cv.size(); i++)
    if (m_split_cv[i].volume > 0) // skip inactive cells
    {
      const auto &cv = m_split_cv[i];
      const size_t parent_dof = m_dof_mapping[i];
      auto &parent_cv = m_cv_data[parent_dof];
      parent_cv.type = cv.type;

      std::cout << "cv = " << i << " parent = " << parent_dof << " ";
      if (cv.type == ControlVolumeType::cell)
        std::cout << "cell" << std::endl;
      else
        std::cout << "face" << std::endl;

      const double volume_fraction = cv.volume / parent_cv.volume;
      parent_cv.aperture += cv.aperture * volume_fraction;
      parent_cv.center += cv.center * volume_fraction;
      parent_cv.porosity += cv.porosity * volume_fraction;
      parent_cv.permeability = cv.permeability; // assume they are the same
      parent_cv.custom = cv.custom;             // assume they are the same
    }
}

void DiscretizationEDFM::build_connection_data_()
{
  create_connections_();
  // we only need to build F-M connections for both dfms and edfms
  // since the frac dofs do not get merged
  for (auto & con : m_con_data)
    if (con.type == ConnectionType::matrix_fracture)
    {
      if (find_edfm_elements_(con).empty())
        build_matrix_dfm_(con);
      else build_matrix_edfm_(con);
    }
}

void DiscretizationEDFM::build_matrix_dfm_(ConnectionData &con)
{
  assert(con.elements.size() == 2);
  // std::cout << con.elements[0] << " " << con.elements[1] << std::endl;
  const auto & cell = m_cv_data[con.elements[1]];
  const auto & frac = m_cv_data[con.elements[0]];
  assert ( cell.type == ControlVolumeType::cell );
  assert ( frac.type == ControlVolumeType::face );
  DiscretizationDFM::build_matrix_fracture(con, frac, cell);
}

void DiscretizationEDFM::build_matrix_edfm_(ConnectionData &con)
{
  assert(con.elements.size() == 2);
  // std::cout << con.elements[0] << " " << con.elements[1] << std::endl;
  const auto & cell = m_cv_data[con.elements[1]];
  const auto & frac = m_cv_data[con.elements[0]];
  assert ( cell.type == ControlVolumeType::cell );
  assert ( frac.type == ControlVolumeType::face );
  const double k_d = con.normal * (cell.permeability * con.normal);  // directional permeability
  const double T = k_d * con.area / con.distances[1]; // average distance from cell to fracture
  con.coefficients = {-T, T};
}

void DiscretizationEDFM::identify_edfm_faces_()
{
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  {
    if (m_edfm_markers.find(face->marker()) != m_edfm_markers.end())
      m_edfm_faces.insert(m_split_dofs.face_dof(face->index()));
  }
}

std::vector<size_t> DiscretizationEDFM::find_edfm_elements_(const ConnectionData & con)
{
  std::vector<size_t> result;
  for (std::size_t i=0; i<con.elements.size(); ++i)
    if (m_edfm_faces.find(con.elements[i]) != m_edfm_faces.end())
      result.push_back(i);
  return result;
}

void DiscretizationEDFM::create_connections_()
{
  hash_algorithms::ConnectionMap<ConnectionData> con_map;
  for (const auto &con : m_split_con)
  {
    const size_t dof1 = m_dof_mapping[con.elements[0]];
    const size_t dof2 = m_dof_mapping[con.elements[1]];
    if ((con.type != ConnectionType::matrix_matrix) && (dof1 == dof2))
    {
      if (!con_map.contains(dof1, dof2)) con_map.insert(dof1, dof2);
      auto &new_con = con_map.get_data(dof1, dof2);
      new_con = con;
      new_con.elements = {dof1, dof2};

      if (con.type == ConnectionType::matrix_fracture)
      {
        const auto & cell = m_split_cv[con.elements[1]];
        const auto & parent_cell = m_cv_data[dof2];
        const auto & face = m_split_cv[con.elements[0]];
        new_con.area += 2 * con.area;
        new_con.normal = con.normal;
        const double dist = (cell.center - face.center).dot( con.normal );
        const double cell_volume_ratio = cell.volume / parent_cell.volume;
        new_con.distances[0] += 0.0;  // from fracture to connection
        new_con.distances[1] += std::fabs(dist) * cell_volume_ratio;  // from cell to connection
      }
    }
  }

  m_con_data = std::move(con_map.get_data());
}

}  // end namespace discretization
