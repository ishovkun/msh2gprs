#include "DiscretizationEDFM.hpp"
#include "DiscretizationDFM.hpp"
#include "DiscretizationTPFA.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;
using mesh::Face;

DiscretizationEDFM::
DiscretizationEDFM(const DoFNumbering & split_dof_numbering,
                   const DoFNumbering & combined_dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data,
                   const std::vector<int> & edfm_markers,
                   const EDFMMethod method)
    : m_split_dofs(split_dof_numbering),
      m_edfm_markers( edfm_markers.begin(), edfm_markers.end() ),
      m_method(method),
      DiscretizationBase(combined_dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationEDFM::build()
{
  std::cout << "build discr" << std::endl;
  DiscretizationDFM discr_dfm(m_split_dofs, m_data, m_split_cv, m_split_con);
  discr_dfm.build();
  if (m_method != EDFMMethod::compartmental)
  {
    identify_edfm_faces_();
    build_control_volume_data_();
    build_connection_data_();
    if (m_method == EDFMMethod::projection)
      build_pedfm_();
    m_con_data = std::move(m_con_map.get_data());
  }
  else
  {
    m_cv_data = std::move(m_split_cv);
    m_con_data = std::move(m_split_con);
  }
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

      // std::cout << "cv = " << i << " parent = " << parent_dof << " ";
      // if (cv.type == ControlVolumeType::cell)
      //   std::cout << "cell" << std::endl;
      // else
      //   std::cout << "face" << std::endl;

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
  std::vector<size_t> m_m_rebuild = create_connections_();
  // we only need to build F-M connections for both dfms and edfms
  // since the frac dofs do not get merged
  std::cout << "rebuild" << std::endl;
  std::cout << "m_con_map.size() = " << m_con_map.size() << std::endl;
  for (const size_t icon : m_m_rebuild)
  {
    auto &con = m_con_map.get_data(icon);
    std::cout << "con = :" ; for (auto i : con.elements) std::cout << i << " ";
    std::cout << std::endl;
    con.center /= con.area;
    DiscretizationTPFA::build_mo(con, m_cv_data[con.elements[0]], m_cv_data[con.elements[1]]);
  }

  std::cout << "new" << std::endl;
  for (auto con = m_con_map.begin(); con != m_con_map.end(); ++con)
    if (con->type == ConnectionType::matrix_fracture)
    {
      if (find_edfm_elements_(*con).empty())
        build_matrix_dfm_(*con);
      else build_matrix_edfm_(*con);
    }
}

void DiscretizationEDFM::build_matrix_dfm_(ConnectionData &con)
{
  assert(con.elements.size() == 2);
  const auto & cell = m_cv_data[con.elements[1]];
  const auto & frac = m_cv_data[con.elements[0]];
  assert ( cell.type == ControlVolumeType::cell );
  assert ( frac.type == ControlVolumeType::face );
  DiscretizationDFM::build_matrix_fracture(con, frac, cell);
}

void DiscretizationEDFM::build_matrix_edfm_(ConnectionData &con)
{
  assert(con.elements.size() == 2);
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

std::vector<size_t> DiscretizationEDFM::create_connections_()
{
  std::vector<size_t> rebuild_connections;
  for (const auto &con : m_split_con)
  {
    const size_t dof1 = m_dof_mapping[con.elements[0]];
    const size_t dof2 = m_dof_mapping[con.elements[1]];
    if (dof1 != dof2)
    {
      if (!m_con_map.contains(dof1, dof2))
      {
        std::cout << "insert " << dof1 << " " << dof2 << std::endl;
        m_con_map.insert(dof1, dof2);
      }

      auto &new_con = m_con_map.get_data(dof1, dof2);
      new_con.elements = {dof1, dof2};

      if ((con.type == ConnectionType::fracture_fracture))
      {
        new_con = con;
      }
      else if (con.type == ConnectionType::matrix_fracture)
      {
        const auto &cell = m_split_cv[con.elements[1]];
        const auto &parent_cell = m_cv_data[dof2];
        const auto &face = m_split_cv[con.elements[0]];
        new_con.area += 2 * con.area;
        new_con.normal = con.normal;
        const double dist = (cell.center - face.center).dot(con.normal);
        const double cell_volume_ratio = cell.volume / parent_cell.volume;
        new_con.distances.resize(2);
        new_con.distances[0] += 0.0; // from fracture to connection
        new_con.distances[1] += std::fabs(dist) * cell_volume_ratio; // from cell to connection
      }
      else // if (con.type == matrix_matrix)
      {
        bool copy_connection = true;
        for (const size_t cv : con.elements)
          if (m_grid.cell(m_split_cv[cv].master).ultimate_parent() != m_grid.cell(m_split_cv[cv].master))
            copy_connection = false;
        if (copy_connection)
        {
          new_con = con;
          new_con.elements = {dof1, dof2};
        }
        else
        {
          rebuild_connections.push_back(m_con_map.index(dof1, dof2));
          new_con.area += con.area;
          new_con.normal = con.normal;
          new_con.center += con.center * con.area;
        }
      }  // end M-M
    }    // end insertion
  }      //end con loop
  return rebuild_connections;
}

void DiscretizationEDFM::build_pedfm_()
{
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  {
    if (m_edfm_markers.find(face->marker()) != m_edfm_markers.end())
    {
      const auto & frac_cell = face->neighbors()[0]->ultimate_parent();
      std::cout << "fine cell = " << face->neighbors()[0]->index() << std::endl;
      std::cout << "frac parent cell = " << frac_cell.index() << std::endl;
      assert(!pedfm_select_faces_(*face).empty());
      for (const mesh::Face* face2 : pedfm_select_faces_(*face))
      {
        const size_t iother_cell = pedfm_find_other_cell_(*face, *face2);
        std::cout << "other parent cell = " << iother_cell << std::endl;
        // have connection between two parent cells
        assert( m_con_map.contains( m_dofs.cell_dof(frac_cell.index()), m_dofs.cell_dof(iother_cell) ) );
        exit(0);
      }
    }

  }
}

std::vector<const mesh::Face*> DiscretizationEDFM::pedfm_select_faces_(const mesh::Face & frac_face) const
{
  //  select largest cell neighbor
  const auto & neighbors = frac_face.neighbors();
  const mesh::Cell* smallest_neighbor;
  if (neighbors[0]->volume() > neighbors[1]->volume())
    smallest_neighbor = neighbors[1];
  else smallest_neighbor = neighbors[0];

  std::vector<const Face*> result;
  std::cout << "smallest_neighbor->faces().size() = " << smallest_neighbor->faces().size() << std::endl;
  int i = 0;
  for (auto face : smallest_neighbor->faces())
  {
    std::cout << i << std::endl;
    i++;
    if (*face == frac_face)                                      // same face
    {
      std::cout << "skip this" << std::endl;
      continue;
    }

    if (face->neighbors().size() < 2)
    {
      std::cout << "skip boundary" << std::endl;
      continue; // domain boundary
    }

    if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-10)  // face ⟂ frac
    {
      std::cout << face->normal() << " | " << frac_face.normal() << std::endl;
      std::cout << "skip perpendicular" << std::endl;
      continue;
    }


    if (m_dofs.is_active_face(face->marker()))
    {
      std::cout << "skip fracs" << std::endl;
      continue; // skip frac faces
    }

    std::cout << "pushing" << std::endl;
    result.push_back(&*face);
  }
  return result;
}

size_t DiscretizationEDFM::pedfm_find_other_cell_(const Face & frac, const Face & other) const
{
  const auto & frac_cell = frac.neighbors()[0]->ultimate_parent();
  for (const auto & cell : other.neighbors())
  {
    const auto & parent = cell->ultimate_parent();
    if (frac_cell != parent)
      return parent.index();
  }
  return frac_cell.index();
}

}  // end namespace discretization
