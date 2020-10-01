#include "DiscretizationPEDFM.hpp"
#include "DiscretizationDFM.hpp"
#include "angem/Collisions.hpp"
#include "Logger.hpp"

namespace discretization {

DiscretizationPEDFM::DiscretizationPEDFM(const DoFNumbering & split_dof_numbering,
                                         const DoFNumbering & combined_dof_numbering,
                                         gprs_data::SimData & data,
                                         std::vector<ControlVolumeData> & cv_data,
                                         std::vector<ConnectionData> & connection_data,
                                         const gprs_data::EmbeddedFractureManager & edfm_mgr)
    : DiscretizationEDFM(split_dof_numbering, combined_dof_numbering,
                         data, cv_data, connection_data,
                         edfm_mgr.get_face_markers()),
      _edfm_mgr(edfm_mgr)
{}

void DiscretizationPEDFM::build()
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
  build_connection_data_();  // regular edfm
  build_non_neighboring_connections_();
  finalize_connections_();
}

void DiscretizationPEDFM::build_non_neighboring_connections_()
{
  // find candidate faces
  for (auto frac_face = m_grid.begin_active_faces(); frac_face != m_grid.end_active_faces(); ++frac_face)
    if (_edfm_mgr.is_fracture(frac_face->marker()))
    {
      const auto isolating_faces = select_faces_(*frac_face);
      const size_t host_cell = host_cell_index_(*frac_face);
      for (const auto * isolating_face : isolating_faces)
        if (non_branching_face_(*frac_face, *isolating_face))
        {
          const size_t neighbor_cell = other_cell_(*frac_face, *isolating_face);
          if (neighbor_cell != host_cell)
          {
            const size_t cell_dof1 = m_dofs.cell_dof(host_cell);
            const size_t cell_dof2 = m_dofs.cell_dof(neighbor_cell);
            const size_t frac_dof = m_dofs.face_dof(frac_face->index());
            // if already cleared, then skip
            if (!_cleared_connections.contains(cell_dof1, cell_dof2))
            {
              const auto projection = project_(_edfm_mgr.fracture_shape(frac_face->marker()),
                                               isolating_face->ultimate_parent().polygon());
              if (projection.first)
              {
                const double projection_area = projection.second->area();
                build_matrix_matrix_(cell_dof1, cell_dof2, projection_area);
                build_fracture_matrix_(m_dofs.face_dof(frac_face->index()), cell_dof2,
                                       projection_area, projection.second->normal());
              }
            }
          }
        }
    }
}

std::list<const mesh::Face*> DiscretizationPEDFM::select_faces_(const mesh::Face & frac_face) const
{
  const mesh::Cell* smallest_neighbor = smallest_neighbor_(frac_face);
  const auto & par = smallest_neighbor->ultimate_parent();

  std::list<const mesh::Face*> result;
  for (auto face : smallest_neighbor->faces())
  {
    if (*face == frac_face) continue; // same face
    if (face->neighbors().size() < 2) continue;  // skip boundary
    if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-6) continue;  // face ⟂ frac
    if (m_dofs.is_active_face(face->index())) continue; // skip frac faces

    result.push_back(face);
  }
  return result;
}

bool DiscretizationPEDFM::non_branching_face_(const mesh::Face & frac_face, const mesh::Face & face) const
{
  // 1. if face does not cross fracture (doesn't have common vertices), then it's not a branch
  if (!have_common_vertices_(frac_face, face))
    return true;
  // 2. if face crosses the fracture and if smaller neighbor is on opposite side of the frac, then
  // it's not a branch
  // 3. otherwise it's a branch
  const auto poly_frac = frac_face.polygon();
  const auto frac_plane = poly_frac.plane();
  const auto cfrac = poly_frac.center();
  const bool side = frac_plane.above(face.center());
  for (const auto * neighbor : face.neighbors())
  {
    const auto neighbor_faces = neighbor->faces();
    // check neighbor cell (that does not contain frac_face)
    if ( std::find(neighbor_faces.begin(), neighbor_faces.end(), &frac_face) == neighbor_faces.end() )
    {
      // a little generic and hard to understand
      // neighbor can be a small sub-cell of a larger host cell
      // we check all sub-cells of the host neighbor cell and divide in regions above and below frac
      // if smaller region of the neighbor cell is on the opposite part of fracture
      // then return true
      if (smaller_cut_part_above_(neighbor, frac_plane) != side)
        return true;
      else return false;
    }
  }

  // this sometimes happens in boundary cells, no idea why
  logging::warning() << "problem in cell " << host_cell_index_(frac_face) << std::endl;
  logging::warning() << "should not be here, it's a bug!" << std::endl;
  // std::cout << "face->neighbors().size() = " << face.neighbors().size() << std::endl;
  // for (const auto * neighbor : face.neighbors())
  //   std::cout << neighbor->index() << " (" << neighbor->ultimate_parent().index() << std::endl;
  // std::cout << "frac? " << _edfm_mgr.is_fracture(face.marker()) << std::endl;
  // std::cout << "marker = " << face.marker() << std::endl;
  // std::cout << "vertices" << std::endl;
  // for (auto v : face.vertices())
  //   std::cout << v << std::endl;
  // // if (face.neighbors().size() == 1)  // skip domain boundaries
  // //   return false;
  // // else
  // throw std::runtime_error("should not be here");
  return false;
}

const mesh::Cell* DiscretizationPEDFM::smallest_neighbor_(const mesh::Face & face) const
{
  const auto & neighbors = face.neighbors();
  if (neighbors.size() < 2) throw std::invalid_argument("face has only one neighbor cell");
  if (neighbors[0]->volume() > neighbors[1]->volume())
    return neighbors[1];
  else return neighbors[0];
}

bool DiscretizationPEDFM::smaller_cut_part_above_(const mesh::Cell* cell,
                                                  const angem::Plane<double> &frac_plane ) const
{
  double vol_above = 0, vol_below = 0;
  for (const auto * subcell : cell->ultimate_parent().ultimate_children())
  {
    if (frac_plane.above(subcell->center())) vol_above += subcell->volume();
    else                                     vol_below += subcell->volume();
  }

  if (vol_above < vol_below) return true;
  else return false;
}

bool DiscretizationPEDFM::have_common_vertices_(const mesh::Face & face1, const mesh::Face & face2) const
{
  auto verts1 = face1.vertices();
  auto verts2 = face2.vertices();
  std::sort(verts1.begin(), verts1.end());
  std::sort(verts2.begin(), verts2.end());
  std::vector<size_t> common_verts;
  std::set_intersection(verts1.begin(), verts1.end(),
                        verts2.begin(), verts2.end(),
                        std::back_inserter(common_verts));
  return !common_verts.empty();
}

size_t DiscretizationPEDFM::other_cell_(const mesh::Face & frac_face,
                                        const mesh::Face & isolating_face) const
{
  const auto & frac_cell = frac_face.neighbors()[0]->ultimate_parent();
  for (const auto & cell : isolating_face.neighbors())
  {
    const auto & parent = cell->ultimate_parent();
    if (frac_cell != parent)
      return parent.index();
  }
  return frac_cell.index();
}

void DiscretizationPEDFM::build_matrix_matrix_(size_t dof1, size_t dof2, double projection_area)
{
  auto & con = m_con_map.get_data(dof1, dof2);
  const auto & cell1 = m_cv_data[dof1];
  const auto & cell2 = m_cv_data[dof2];

  const auto &c1 = cell1.center;
  const auto &c2 = cell2.center;
  const auto &cf = con.center;
  const auto &n = con.normal;

  const double t = (cf - c1).dot(n) / (c2 - c1).dot(n);
  const auto cp = c1 + t * (c2 - c1);
  const auto &K1 = cell1.permeability;
  const auto &K2 = cell2.permeability;
  const double Kp1 = (K1 * (c1 - cp).normalize()).norm();
  const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
  const double subtracted_area = projection_area;
  const double dT1 = subtracted_area * Kp1 / (c1 - cp).norm();
  const double dT2 = subtracted_area * Kp2 / (c2 - cp).norm();
  const double dT = (dT1 + dT2 > 0) ? dT1 * dT2 / (dT1 + dT2) : 0.0;
  const double T_new = std::max(std::fabs(con.coefficients[0]) - dT, 0.0);
  con.coefficients = {-T_new, T_new};

  // compute initial trans to decide whether to kill the connection
  const double T1_init = con.area * Kp1 / c1.distance(cp);
  const double T2_init = con.area * Kp2 / c2.distance(cp);
  const double T_init = (T1_init + T2_init > 0) ? T1_init * T2_init / (T1_init + T2_init) : 0.0;
  std::cout << T_new << "/" << T_init;
  if (T_new < 1e-6 * T_init)
  {
    std::cout << " clear";
    _cleared_connections.insert(dof1, dof2);
  }
  std::cout << std::endl;
}

void DiscretizationPEDFM::build_fracture_matrix_(size_t frac_dof, size_t cell_dof,
                                                 double projection_area,
                                                 const angem::Point<3,double> & normal)
{
  const size_t con_idx = m_con_map.insert(frac_dof, cell_dof);
  auto & con = m_con_map.get_data(con_idx);
  con.normal = normal;
  con.area = projection_area;
  con.elements = {frac_dof, cell_dof};

  const auto & frac = m_cv_data[frac_dof];
  const auto & cell = m_cv_data[cell_dof];

  const auto &c1 = frac.center;
  const auto &c2 = cell.center;
  const auto &cf = con.center;
  const auto &n = con.normal;
  // projection point
  const double t =  (cf - c1).dot(n) / (c2 - c1).dot(n);
  const auto cp = c1 + t*(c2 - c1);
  // project permeability
  const double & Kp1 = frac.permeability(0, 0);;
  const auto & K2 = cell.permeability;
  const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
  // cell-face transmissibility
  const double face_area = con.area;
  const double T1 = face_area * Kp1 / (c1 - cp).norm();
  const double T2 = face_area * Kp2 / (c2 - cp).norm();
  const double distance_tot = (c1 - c2).norm();
  // face transmissibility
  double T = 0.0;
  if ( std::isnormal(T1 + T2) )
    T = T1*T2 / ( T1 + T2 );
  con.coefficients = {-T, T};
  con.update_formula = {face_area, distance_tot, Kp2, T };
}

void DiscretizationPEDFM::finalize_connections_()
{
  // copy to output vector while avoiding to copy cleared connections condata
  m_con_data.reserve(m_con_map.size());
  for (auto it = m_con_map.begin(); it != m_con_map.end(); ++it)
    if (!_cleared_connections.contains(it.elements()))
        m_con_data.push_back(*it);
}

std::pair<bool, std::unique_ptr<angem::Polygon<double>>>
DiscretizationPEDFM::project_(const angem::Polygon<double> & poly1, const angem::Polygon<double> & poly2) const
{
  std::pair<bool, std::unique_ptr<angem::Polygon<double>>> result;
  const double tol = 1e-6;
  if (std::fabs(poly1.normal().dot( poly2.normal() )) < tol)
  {
    result.first = false;
    return result;
  }

  const auto projected_frac_vertices = poly2.plane().project_points(poly1.get_points());
  const angem::Polygon<double> proj_poly(projected_frac_vertices);
  std::vector<angem::Point<3,double>> intersection;
  result.first = angem::collision(proj_poly, poly2, intersection, tol);
  if (result.first)
    result.second = std::make_unique<angem::Polygon<double>>(intersection);

  return result;
}

size_t DiscretizationPEDFM::host_cell_index_(const mesh::Face & frac_face) const
{
  return frac_face.neighbors()[0]->ultimate_parent().index();}

}  // end namespace discretization
