#include "DiscretizationTPFA.hpp"
#include "angem/Point.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3,double>;

DiscretizationTPFA::
DiscretizationTPFA(const DoFNumbering & dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data)
    :
    DiscretizationBase(dof_numbering, data, cv_data, connection_data),
    m_method(tpfa_method::mo)
{}


void DiscretizationTPFA::build()
{
  DiscretizationBase::build_cell_data_();

  m_con_data.reserve(m_con_data.size() + m_grid.n_faces());
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  {
    if (face->neighbors().size() == 2)
      if (!m_dofs.is_active_face(face->index()))
      {
        m_con_data.emplace_back();
        auto & con = m_con_data.back();
        const auto cells = face->neighbors();
        const size_t cv1 = m_dofs.cell_dof(cells[0]->index());
        const size_t cv2 = m_dofs.cell_dof(cells[1]->index());
        con.elements = {cv1, cv2};
        con.center = face->center();
        con.normal = face->normal();
        con.area = face->area();
        con.type = ConnectionType::matrix_matrix;
        build_mo(con, m_cv_data[ con.elements[0] ], m_cv_data[ con.elements[1] ]);
      }
  }
}


void DiscretizationTPFA::build_kirill(const mesh::Face & face,
                                      ConnectionData   & data)
{
  // Point nfmf, center0, center1, center_face;
  const auto cells = face.neighbors();
  const size_t cv1 = m_dofs.cell_dof(cells[0]->index());
  const size_t cv2 = m_dofs.cell_dof(cells[1]->index());
  const auto & cell1 = m_cv_data[ cv1 ];
  const auto & cell2 = m_cv_data[ cv2 ];

  const Point center_face = face.center();
  const Point d1 = cell1.center - center_face;
  const Point d2 = cell2.center - center_face;

  const Point face_normal = face.normal();
  const Tensor & K1 = cell1.permeability;
  const Tensor & K2 = cell2.permeability;

  // project permeability on face normal
  const Point K1_n = K1*face_normal;
  const Point K2_n = K2*face_normal;

  // projection of distance along normal
  const double dd1 = fabs(face_normal.dot(d1));
  const double dd2 = fabs(face_normal.dot(d2));

  // projection of conormal along normal
  const double n_K1_n = face_normal.dot(K1_n);
  const double n_K2_n = face_normal.dot(K2_n);

  double T; // transmissibility
  double D; // projection distance
  static const double eps = 1e-8;
  if (dd1 > eps && dd2 > eps)
  {
    T = n_K1_n * n_K2_n / (dd1*n_K2_n + dd2*n_K1_n + 1.0e-35);
    D = 1.0 / (dd1 + dd2 + 1.0e-35);
  }
  else if (dd1 > eps)
  {
    T = n_K1_n / dd1;
    D = 1.0 / dd1;
  }
  else if (dd2 > eps)
  {
    T = n_K2_n / dd2;
    D = 1.0 / dd2;
  }
  else
  {
    T = 0;
    D = 0.0;
  }

  // T /= D;
  T *= face.area();

  // NOTE: We need to save D for the TPFACONNSN keyword

  data.elements.resize(2);
  data.elements[0] = cv1;
  data.elements[1] = cv2;

  data.coefficients.resize(2);
  data.coefficients[0] = -T;
  data.coefficients[1] =  T;
}


void DiscretizationTPFA::build_mo(ConnectionData & con,
                                  const ControlVolumeData & cell1,
                                  const ControlVolumeData & cell2)
{

  // define face projection point
  const Point & c1 = cell1.center;
  const Point & c2 = cell2.center;
  const Point & cf = con.center;
  const Point & n = con.normal;
  // projection point
  const double t =  (cf - c1).dot(n) / (c2 - c1).dot(n);
  const Point cp = c1 + t*(c2 - c1);
  // project permeability
  const Tensor & K1 = cell1.permeability;
  const Tensor & K2 = cell2.permeability;
  // directional permeability
  const double Kp1 = (K1 * (c1 - cp).normalize()).norm();
  const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
  // cell-face transmissibility
  const double face_area = con.area;
  const double T1 = face_area * Kp1 / (c1 - cp).norm();
  const double T2 = face_area * Kp2 / (c2 - cp).norm();
  // face transmissibility
  const double T = T1*T2 / ( T1 + T2 );
  con.coefficients = {-T, T};
}

}
