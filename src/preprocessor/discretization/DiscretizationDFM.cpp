#include "DiscretizationDFM.hpp"
#include "angem/Tensor2.hpp"
#include <cmath>  // std::isnan
#include <numeric>  // std::accumulate

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3,double>;

DiscretizationDFM::
DiscretizationDFM(const mesh::Mesh                              & grid,
                  const std::set<int>                           & dfm_markers,
                  const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                  const std::vector<std::vector<double>>        & props,
                  const std::vector<std::string>                & keys,
                  const size_t                                    shift_matrix,
                  const size_t                                    shift_dfm)
: DiscretizationBase(grid, dfm_markers, props, keys),
  dfm_faces(dfm_faces), shift_matrix(shift_matrix),
  shift_dfm(shift_dfm)
{
  // check that ranges don't intersect
  assert( (shift_matrix + grid.n_cells() <= shift_dfm) ||
          (shift_matrix >= shift_dfm + dfm_faces.size()) );
}


void DiscretizationDFM::build()
{
  std::cout << "build cv data" << std::endl;
  build_cell_data();

  // build connection lists (no data)
  build_fracture_matrix_connections();

  std::cout << "build F-M discretization" << std::endl;
  for (auto & con : con_data)
    build_matrix_fracture(con);

  // build fracture-fracture transes
  std::cout << "build F-F discretization" << std::endl;
  build_fracture_fracture_connections();
}


void DiscretizationDFM::build_cell_data()
{
  // could be less though, just to be save
  cell_to_cv.resize(grid.n_cells());
  cv_data.resize(dfm_faces.size() * 3);
  std::unordered_set<size_t> bounding_cells;
  apertures.resize(dfm_faces.size());

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 2)
    {
      if (is_fracture(face.marker()))
      {
        const auto it = dfm_faces.find(face.index());
        assert(it != dfm_faces.end());
        const auto & face_props = it->second;

        auto & data = cv_data[face_props.cv_index];
        data.type = ControlVolumeType::face;
        std::cout << "dfm face " << face_props.cv_index << std::endl;
        data.master = face.index();
        data.center = face.center();
        data.volume = face_props.aperture * face.area();
        data.permeability = Tensor::make_unit_tensor();
        data.permeability *= (face_props.conductivity / face_props.aperture);
        data.porosity = 1.0;
        apertures[face_props.cv_index] = face_props.aperture;

        // take custom properties as weighted average from bounding cells
        const auto cells = face.neighbors();
        const auto cell1 = grid.create_const_cell_iterator(cells[0]);
        const auto cell2 = grid.create_const_cell_iterator(cells[1]);
        data.custom.resize(custom_keys.size());
        for (size_t j = 0; j < custom_keys.size(); ++j)
          data.custom[j] +=
              (props[cell1.index()][custom_keys[j]] * cell1.volume() +
               props[cell2.index()][custom_keys[j]] * cell2.volume()) /
              ( cell1.volume() + cell2.volume() ) ;

        bounding_cells.insert(cells[0]);
        bounding_cells.insert(cells[1]);
      }
    }

  size_t cv = dfm_faces.size();
  for (const std::size_t i : bounding_cells)
  {
    std::cout << i << std::endl;
    const auto cell = grid.create_const_cell_iterator(i);
    auto & data = cv_data[cv];
    // auto & data = cv_data.emplace_back();
    cell_to_cv[i] = cv;

    data.type = ControlVolumeType::cell;
    data.master = i;
    data.center = cell.center();
    data.volume = cell.volume();
    data.permeability = get_permeability(i);

    data.custom.resize(custom_keys.size());
    for (size_t j = 0; j < custom_keys.size(); ++j)
      data.custom[j] = props[i][custom_keys[j]];

    cv++;
  }

  cv_data.resize(cv);
}


hash_algorithms::ConnectionMap<std::vector<size_t>>
DiscretizationDFM::map_edge_to_faces()
{
  hash_algorithms::ConnectionMap<std::vector<size_t>> edge_face_connections;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 2)
      if (is_fracture(face.marker()))
      {
        const auto it = dfm_faces.find(face.index());
        assert(it != dfm_faces.end());
        const auto & props = it->second;

        for (const auto & edge : face.edges())
        {
          size_t index;
          if (edge_face_connections.has(edge.first, edge.second))
            index = edge_face_connections.index(edge.first, edge.second);
          else
            index = edge_face_connections.insert(edge.first, edge.second);

          auto & cvs = edge_face_connections.get_data(index);

          // get frac properties
          cvs.push_back(props.cv_index);
        }
      }
  return edge_face_connections;
}


void DiscretizationDFM::build_fracture_fracture(ConnectionData & con)
{
  std::cout << "not implemented" << std::endl;
  abort();
}


void DiscretizationDFM::build_matrix_fracture(ConnectionData & con)
{
  // cause I built them that way
  const auto & cv_frac = cv_data[con.elements[0]];
  const auto & cv_cell = cv_data[con.elements[1]];
  assert(cv_frac.type == ControlVolumeType::face);
  assert(cv_cell.type == ControlVolumeType::cell);
  assert(con.type == ConnectionType::matrix_fracture);

  // project cell permeability
  const auto f = cv_cell.center - con.center;
  const double K_cell = (cv_cell.permeability * (f/f.norm())).norm();

  // frac perm is just conductivilty / aperture
  const double K_frac = cv_frac.permeability(0, 0);

  const double T_cell = con.area * K_cell / f.norm();
  const double T_face = con.area * K_frac;

  // connection transmissibility
  double T = 0;
  if ( !std::isnan(1. / (T_cell + T_face) ) )
    T = T_cell*T_face / (T_cell + T_face);

  con.coefficients.resize(2);
  con.coefficients[0] = -T;
  con.coefficients[1] =  T;
}


void DiscretizationDFM::build_fracture_matrix_connections()
{
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    if (is_fracture(face.marker()))
    {
      const auto & neighbors = face.neighbors();
      // find frac cv index
      std::cout << face.index() << std::endl;
      const auto it = dfm_faces.find(face.index());
      assert(it != dfm_faces.end());
      const std::size_t cv_frac = it->second.cv_index;

      const auto & face_props = it->second;
      // connection fracture-cell1
      {
        con_data.emplace_back();
        auto &con = con_data.back();
        con.elements.push_back(cv_frac);
        con.elements.push_back(cell_to_cv[neighbors[0]]);
        con.type = ConnectionType::matrix_fracture;
        con.area = face.area();
        con.normal = face.normal();
        con.center = face.center();
      }

      //  connection fracture-cell2
      {
        con_data.emplace_back();
        auto &con = con_data.back();
        con.type = ConnectionType::matrix_fracture;
        con.elements.push_back(cv_frac);
        con.elements.push_back(cell_to_cv[neighbors[1]]);
        con.area = face.area();
        con.normal = face.normal();
        con.center = face.center();
      }
    }
  }
}


void DiscretizationDFM::build_fracture_fracture_connections()
{
  // build fracture-fracture-connections
  // map edges to dfm faces
  auto edge_face_connections = map_edge_to_faces();
  // we won't need all of those since there are edges connected to
  // only one face
  con_data.reserve(2*dfm_faces.size() + edge_face_connections.size());
  for (auto edge = edge_face_connections.begin();
       edge != edge_face_connections.end(); ++edge)
  {
    const std::vector<size_t> & face_cvs = *edge;
    if (face_cvs.size() > 1)
    {
      const auto edge_vertices = edge.elements();
      const Point e1 = grid.vertex_coordinates(edge_vertices.first);
      const Point e2 = grid.vertex_coordinates(edge_vertices.second);
      const Point edge_center = 0.5 * (e1 + e2) ;
      const Point de = e2 - e1;
      const double edge_length = de.norm();

      // compute average (by number) projection onto the edge
      Point cv_projection = edge_center;
      for (std::size_t i = 0; i < face_cvs.size(); ++i)
      {
        const double t = de.dot(cv_data[face_cvs[i]].center - edge_center) / edge_length;
        cv_projection += t * de / face_cvs.size();
      }

      // compute parts of transmissibility
      std::vector<double> transmissibility_part(face_cvs.size());
      for (std::size_t i = 0; i < face_cvs.size(); ++i)
      {
        // aperture * edge length
        const double area = apertures[face_cvs[i]] * edge_length;
        const double dist_to_edge = (cv_data[face_cvs[i]].center - edge_center).norm();
        const double perm = cv_data[face_cvs[i]].permeability(0, 0);
        transmissibility_part[i] = area * perm / dist_to_edge;
      }

      const double t_sum = std::accumulate(transmissibility_part.begin(),
                                           transmissibility_part.end(), 0);
      for (std::size_t i = 0; i < face_cvs.size(); ++i)
        for (std::size_t j = i+1; j < face_cvs.size(); ++j)
        {
          auto & con = con_data.emplace_back();
          con.type = ConnectionType::fracture_fracture;
          con.elements = {face_cvs[i], face_cvs[j]};
          con.coefficients.resize(2);
          con.coefficients[0] = transmissibility_part[i] *
                                transmissibility_part[j] / t_sum;
          con.coefficients[1] = -con.coefficients[0];
        }
    }
  }

}

} // end namespace
