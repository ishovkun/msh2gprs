#include "GridGeneratorINSIM.hpp"
#include "Well.hpp"                                       // provides Well
#include "gmsh_interface/GmshInterface.hpp"
#ifdef WITH_GMSH
#include <gmsh.h>
#endif
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter

namespace gprs_data {

GridGeneratorINSIM::GridGeneratorINSIM(INSIMMeshConfig const & config, std::vector<WellConfig> const & wells)
    : _config(config), _wells(wells)
{
  if ( _wells.size() < 2 )
    throw std::invalid_argument("Cannot build INSIM grid with less than two wells");

  for (auto const & well : _wells) {
    assert( well.coordinates.size() == 1 && "Only simple wells are supported for now" );
    setup_simple_well_(well);
  }

  generate_bounding_box_();
}

void GridGeneratorINSIM::generate_bounding_box_()
{
  double const margin = find_characteristic_length_();
  double const upper = std::numeric_limits<double>::max();
  double const lower = std::numeric_limits<double>::lowest();
  angem::Point<3,double> bbox_min = {upper, upper, upper};
  angem::Point<3,double> bbox_max = {lower, lower, lower};
  for (auto const & v : _bounds)
    for (size_t i = 0; i < 3; ++i) {
      bbox_min[i] = std::min( bbox_min[i], v[i] );
      bbox_max[i] = std::max( bbox_max[i], v[i] );
    }

  bbox_min -= margin;
  bbox_max += margin;

  auto const delta = bbox_max - bbox_min;

  std::vector<angem::Point<3,double>> verts(8);  // hex has 8 vertices
  std::fill( verts.begin(), verts.begin() + 4, bbox_min );

  // look at the schematics in angem/Hexahedron.hpp to understand the order
  verts[1][0] += delta[0];
  verts[2][0] += delta[0];
  verts[2][1] += delta[1];
  verts[3][1] += delta[1];

  // points 4-7 are the same but shifted by z
  std::copy( verts.begin(), verts.begin() + 4, verts.begin() + 4 );
  std::for_each( verts.begin() + 4, verts.end(), [&delta](auto & v) { v[2] += delta[2];} );
  // range array
  std::vector<std::size_t> indices(verts.size());
  std::iota( indices.begin(), indices.end(), 0 );
  // make hexahedron
  _bbox = std::make_unique<angem::Hexahedron<double>>(verts, indices);
}

double GridGeneratorINSIM::find_characteristic_length_() const
{
  // compute the center of all well nodes
  angem::Point<3,double> const c = angem::compute_center_mass( _bounds );

  // compute average distance between c and vertices
  double dist = 0.f;
  for (auto const & v : _bounds)
    dist += v.distance(c);

  return dist / _bounds.size();
}


GridGeneratorINSIM::operator mesh::Mesh() const
{
  mesh::Mesh grid;

#ifdef WITH_GMSH
  GmshInterface::initialize_gmsh(/*verbose = */ true);
  GmshInterface::build_triangulation_embedded_points(*_bbox, _embedded_pts, grid);
GmshInterface::finalize_gmsh();
#endif
  return grid;
}

void GridGeneratorINSIM::setup_simple_well_(WellConfig const & conf)
{
  assert( _config.minimum_thickness > 0 );
  auto v1 = conf.coordinates[0];
  v1[2] += 0.5 * _config.minimum_thickness;  // shift up in z direction
  _bounds.push_back( v1 );
  auto v2 = v1;
  v2[2] -= 0.5 * _config.minimum_thickness;  // shift down in z direction
  _bounds.push_back( v2 );
  size_t const vertex_idx = _embedded_pts.size();
  // _well_vertices.push_back({vertex_idx});
  _embedded_pts.push_back(conf.coordinates[0]);
}

}  // end namespace gprs_data
