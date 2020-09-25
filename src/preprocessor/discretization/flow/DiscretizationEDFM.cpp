#include "DiscretizationEDFM.hpp"
#include "DiscretizationDFM.hpp"
#include "DiscretizationTPFA.hpp"
#include "angem/Projections.hpp"  // provides angem::project

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3, double>;
using mesh::Face;

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
  if (m_edfm_faces.empty())
  {
    m_cv_data = std::move(m_split_cv);
    m_con_data = std::move(m_split_con);
    return;
  }

  build_control_volume_data_();
  build_connection_data_();

  // if (m_method == EDFMMethod::projection)
  //   build_pedfm_();

  // copy to condata
  m_con_data.reserve(m_con_map.size());
  for (auto it = m_con_map.begin(); it != m_con_map.end(); ++it)
    m_con_data.push_back(*it);
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
  {
    const auto &cv = m_split_cv[i];
    size_t parent_dof;
    if (cv.type == ControlVolumeType::cell)
    {
      parent_dof = m_dofs.cell_dof(cv.master);
      m_cv_data[parent_dof].master = m_grid.cell(cv.master).ultimate_parent().index();
    }
    else // if (cv.type == ControlVolumeType::face)
    {
      parent_dof = m_dofs.face_dof(cv.master);
      m_cv_data[parent_dof].master = cv.master;
    }

    m_dof_mapping[i] = parent_dof;
    m_cv_data[parent_dof].volume += cv.volume;
  }

  for (size_t i = 0; i < m_split_cv.size(); i++)
    {
      const auto &cv = m_split_cv[i];
      const size_t parent_dof = m_dof_mapping[i];
      auto &parent_cv = m_cv_data[parent_dof];
      parent_cv.type = cv.type;

      const double volume_fraction = cv.volume / parent_cv.volume;
      parent_cv.aperture += cv.aperture * volume_fraction;
      parent_cv.center += cv.center * volume_fraction;
      parent_cv.porosity += cv.porosity * volume_fraction;
      parent_cv.permeability = cv.permeability; // assume they are the same
      parent_cv.custom = cv.custom;             // assume they are the same
    }

  for(size_t i = 0; i < m_cv_data.size(); i++)
  {
    const auto & cv = m_cv_data[i];
    // std::cout << "cv.index= " << i << std::endl;
    if ( cv.type == ControlVolumeType::cell && cv.volume == 0 )
      throw std::runtime_error("discr edfm zero cell volume " +
                               std::to_string(i) + " (cell " +
                               std::to_string(cv.master) + ")");
    else if ( cv.type == ControlVolumeType::face && cv.volume == 0 )
      throw std::runtime_error("discr edfm zero face volume");
  }

}

void DiscretizationEDFM::build_connection_data_()
{
  std::vector<size_t> m_m_rebuild = create_connections_();
  // we only need to build F-M connections for both dfms and edfms
  // since the frac dofs do not get merged
  for (const size_t icon : m_m_rebuild)
  {
    auto &con = m_con_map.get_data(icon);
    con.center /= con.area;
    assert( con.elements[0] < m_cv_data.size() );
    assert( con.elements[1] < m_cv_data.size() );
    DiscretizationTPFA::build_mo(con, m_cv_data[con.elements[0]], m_cv_data[con.elements[1]]);
  }

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
  const double d = con.distances[1]; // average distance from cell to fracture
  const double k_d = con.normal * (cell.permeability * con.normal);  // directional permeability
  const double k_f = frac.permeability(0, 0);                        // frac permeability
  const double T_m = k_d * con.area / d; // matrix half
  const double T_f = k_f * con.area / d; // matrix half
  //
  //  connection transmissibility
  double T = 0;
  if ( !std::isnan(1. / (T_m + T_f) ) )
    T = T_m*T_f / (T_m + T_f);
  con.coefficients = {-T, T};
  con.update_formula = { con.area, d, k_d, T };
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
  std::unordered_set<size_t> rebuild_connections;
  for (const auto &con : m_split_con)
  {
    assert( con.elements[0] < m_dof_mapping.size() );
    assert( con.elements[1] < m_dof_mapping.size() );
    const size_t dof1 = m_dof_mapping[con.elements[0]];
    const size_t dof2 = m_dof_mapping[con.elements[1]];

    if (dof1 != dof2)
    {
      if (!m_con_map.contains(dof1, dof2))
      {
        // std::cout << "insert " << m_cv_data[dof1].master << " " <<m_cv_data[dof2].master
        //     << " (" << dof1 << " " << dof2 << ") "<< std::endl;
        m_con_map.insert(dof1, dof2);
        // set needed to zero
        auto &new_con = m_con_map.get_data(dof1, dof2);
        new_con.type = con.type;
        new_con.elements = {dof1, dof2};
      }

      auto &new_con = m_con_map.get_data(dof1, dof2);

      if ((con.type == ConnectionType::fracture_fracture))
      {
        // since we don't merge split fracture elements
        new_con = con;
      }
      else if (con.type == ConnectionType::matrix_fracture)
      {
        const auto &cell = m_split_cv[con.elements[1]];
        const auto &parent_cell = m_cv_data[dof2];
        const auto &face = m_split_cv[con.elements[0]];

        if (!find_edfm_elements_(con).empty())  // EDFM-matrix
        {
          new_con.area = 2 * con.area;
          const double dist = (cell.center - face.center).dot(con.normal);
          const double cell_volume_ratio = cell.volume / parent_cell.volume;
          // Calculate the volumetric average
          if(new_con.distances.empty())
            new_con.distances = {0.0, 0.0};
          new_con.distances[1] += std::fabs(dist) * cell_volume_ratio;
        }
        else // DFM-matrix
        {
          new_con.area = con.area;
        }

        // common for EDFM-M and DFM-M
        new_con.normal = con.normal;
        new_con.center = face.center;

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
          rebuild_connections.insert(m_con_map.index(dof1, dof2));
          new_con.area += con.area;
          new_con.normal = con.normal;
          new_con.center += con.center * con.area;
        }
      }  // end M-M
    }    // end insertion
  }      //end con loop
  return std::vector<size_t>(rebuild_connections.begin(), rebuild_connections.end());
}

void DiscretizationEDFM::build_pedfm_()
{
  std::cout << "\n\nstuff:" << std::endl;
  PureConnectionMap cleared_connections;
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  {
    if (m_edfm_markers.find(face->marker()) != m_edfm_markers.end())
    {
      const auto & frac_cell = face->neighbors()[0]->ultimate_parent();
      // std::cout << "searching nighbors for frac cell " << frac_cell.index() << std::endl;

      assert(!pedfm_select_faces_(*face).empty());
      for (const mesh::Face* face2 : pedfm_select_faces_(*face))
      {
        const size_t iother_cell = pedfm_find_other_cell_(*face, *face2);
        // determing all participating dofs
        const size_t cell_dof1 = m_dofs.cell_dof(frac_cell.index());
        const size_t cell_dof2 = m_dofs.cell_dof(iother_cell);
        const size_t face_dof = m_dofs.face_dof(face->index());
        // if already cleared, then skip
        if (cleared_connections.contains(cell_dof1, cell_dof2))
          continue;
        // this happens when several edfm's in a cell
        if (iother_cell == frac_cell.index()) continue;

        // have connection between two parent cells
        assert( m_con_map.contains( cell_dof1, cell_dof2 ) );

        try
        {
          // projection of frac face onto M-M connecting face
          const auto projection_points = angem::project(face->polygon(), face2->polygon(), 1e-6);
          if (projection_points.size() < 3)
            continue; // no real projection

          const size_t ifrac_mat_con = m_con_map.insert(face_dof, cell_dof2);
          auto &frac_mat_con = m_con_map.get_data(ifrac_mat_con);
          frac_mat_con.elements = {face_dof, cell_dof2};
          // note that we invoke it after the insertion insertion killls
          // references
          auto &mat_mat_con = m_con_map.get_data(cell_dof1, cell_dof2);

          const angem::Polygon<double> projection(projection_points);
          frac_mat_con.normal = projection.normal();
          frac_mat_con.area = projection.area();
          const bool kill_connection = build_pedfm_(mat_mat_con, frac_mat_con);
          // if small trans then kill connection
          if (kill_connection) {
        // if (frac_cell.index() == 1482 || iother_cell == 1482)
        //     {
        //       std::cout << "clear connection " << frac_cell.index() << " "
        //                 << iother_cell << "(" << cell_dof1 << ", "<< cell_dof2 << ")" <<std::endl;
        //     }

            cleared_connections.insert(cell_dof1, cell_dof2);
            // exit(0);
          }
        }
        catch(const std::runtime_error & err)
        {
          std::cout << "zero projection" << std::endl;
          // somtimees the projection on the face plane
          // is outside of the M-M face. This usually
          // happens when multiple embedded fractures cross a cell
          continue;
        }

      }
    }
  }

  // clear connections
  // std::cout << "clean con" << std::endl;
  for (auto it = cleared_connections.begin(); it!=cleared_connections.end(); ++it)
  {
    const auto pair_elements = it.elements();
    // std::cout << "clear " <<  pair_elements.first << " " <<  pair_elements.second << std::endl;
    m_con_map.remove(pair_elements.first, pair_elements.second);
  }
  std::cout << "\n\n\n" << std::endl;
}

std::vector<const mesh::Face*> DiscretizationEDFM::pedfm_select_faces_(const mesh::Face & frac_face) const
{
  const auto & neighbors = frac_face.neighbors();

  //  method 1 by volume
  //    select largest cell neighbor
  const mesh::Cell* smallest_neighbor;
  if (neighbors[0]->volume() > neighbors[1]->volume())
    smallest_neighbor = neighbors[1];
  else smallest_neighbor = neighbors[0];
  const auto & par = smallest_neighbor->ultimate_parent();

  std::vector<const Face*> result;
  for (auto face : smallest_neighbor->faces())
  {
    if (*face == frac_face) continue; // same face
    if (face->neighbors().size() < 2) continue;  // skip boundary
    if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-6) continue;  // face ⟂ frac
    if (m_dofs.is_active_face(face->index())) continue; // skip frac faces

    result.push_back(&*face);
  }
  return result;

  // method 2
  // const auto & par = neighbors[0]->ultimate_parent();


  // // select closest and pick closest neighbor
  // std::vector<const Face*> result1;
  // double dist1 = std::numeric_limits<double>::max();
  // for (auto face : neighbors[0]->faces())
  // {
  //   if (*face == frac_face) continue; // same face
  //   if (face->neighbors().size() < 2) continue;  // skip boundary
  //   if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-6) continue;  // face ⟂ frac
  //   if (m_dofs.is_active_face(face->index())) continue; // skip frac faces
  //   dist1 = std::min( dist1, frac_face.center().distance(face->center()) );

  //   result1.push_back(&*face);
  // }

  // std::vector<const Face*> result2;
  // double dist2 = std::numeric_limits<double>::max();
  // for (auto face : neighbors[1]->faces())
  // {
  //   if (*face == frac_face) continue; // same face
  //   if (face->neighbors().size() < 2) continue;  // skip boundary
  //   if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-6) continue;  // face ⟂ frac
  //   if (m_dofs.is_active_face(face->index())) continue; // skip frac faces
  //   dist2 = std::min( dist2, frac_face.center().distance(face->center()) );

  //   result2.push_back(&*face);
  // }

  // if (dist1 < dist2)
  //   return result1;
  // else return result2;
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

bool DiscretizationEDFM::build_pedfm_(ConnectionData & mm_con,
                                      ConnectionData & fm_con)
{
  assert(mm_con.elements.size() == 2);
  assert(fm_con.elements.size() == 2);
  auto pair_cells = find_fracture_cv_and_nonfracture_cv_(mm_con,fm_con);
  assert( pair_cells.first < m_cv_data.size() );
  assert( pair_cells.second < m_cv_data.size() );
  assert( fm_con.elements[0] < m_cv_data.size() );

  const auto & cell1 = m_cv_data[pair_cells.first];
  const auto & cell2 = m_cv_data[pair_cells.second];
  const auto & frac = m_cv_data[fm_con.elements[0]];

  bool kill_connection = false;

  // first modify old M_M connection
  {
    const Point &c1 = cell1.center;
    const Point &c2 = cell2.center;
    const Point &cf = mm_con.center;
    const Point &n = mm_con.normal;
    const double t = (cf - c1).dot(n) / (c2 - c1).dot(n);
    const Point cp = c1 + t * (c2 - c1);
    const Tensor &K1 = cell1.permeability;
    const Tensor &K2 = cell2.permeability;
    const double Kp1 = (K1 * (c1 - cp).normalize()).norm();
    const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
    const double subtracted_area = fm_con.area;
    const double dT1 = subtracted_area * Kp1 / (c1 - cp).norm();
    const double dT2 = subtracted_area * Kp2 / (c2 - cp).norm();
    const double dT = (dT1 + dT2 > 0) ? dT1 * dT2 / (dT1 + dT2) : 0.0;
    const double T_new = std::max(std::fabs(mm_con.coefficients[0]) - dT, 0.0);
    mm_con.coefficients = {-T_new, T_new};

    // compute initial trans to decide whether to kill the connection
    const double T1_init = mm_con.area * Kp1 / c1.distance(cp);
    const double T2_init = mm_con.area * Kp2 / c2.distance(cp);
    const double T_init = (T1_init + T2_init > 0) ? T1_init * T2_init / (T1_init + T2_init) : 0.0;
    if (cell1.master == 1482 || cell2.master == 1482)
    {
      std::cout << "areas: "<< subtracted_area << " " << mm_con.area << std::endl;
      std::cout << "trans: " << T_new << " " << T_init << std::endl;

    }
    if (T_new < 1e-6 * T_init)
    {
      // std::cout << "subtracted = " << subtracted_area << std::endl;
      // std::cout << "mm_con.area = " << mm_con.area << std::endl;
      // std::cout << "ratio = " << new_mm_trans / T_init << std::endl;
      kill_connection = true;
    }

  }

  // now compute non-neighboring F-M trans
  {
    const Point &c1 = frac.center;
    const Point &c2 = cell2.center;
    const Point &cf = fm_con.center;
    const Point &n = fm_con.normal;
    // projection point
    const double t =  (cf - c1).dot(n) / (c2 - c1).dot(n);
    const Point cp = c1 + t*(c2 - c1);
    // project permeability
    const Tensor & K2 = cell2.permeability;
    const double & Kp1 = frac.permeability(0, 0);;
    const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
    // cell-face transmissibility
    const double face_area = fm_con.area;
    const double T1 = face_area * Kp1 / (c1 - cp).norm();
    const double T2 = face_area * Kp2 / (c2 - cp).norm();
    const double distance_tot = (c1 - c2).norm();
    // face transmissibility
    double T = 0.0;
    if ( std::isnormal(T1 + T2) )
      T = T1*T2 / ( T1 + T2 );
    fm_con.coefficients = {-T, T};
    fm_con.update_formula = {face_area, distance_tot, Kp2, T };
  }
  return kill_connection;
}

std::pair<size_t,size_t> DiscretizationEDFM::find_fracture_cv_and_nonfracture_cv_(const ConnectionData & mm_con, const ConnectionData & fm_con) const
{
  if ( fm_con.elements[1] == mm_con.elements[0] )
    return {mm_con.elements[0], mm_con.elements[1]};
  else return {mm_con.elements[1], mm_con.elements[0]};
}

}  // end namespace discretization
