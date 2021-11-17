#include "GridGeneratorINSIM.hpp"
#include "Well.hpp"                                       // provides Well
#include "gmsh_interface/GmshInterface.hpp"
#ifdef WITH_GMSH
#include <gmsh.h>
#endif
// #include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter

namespace gprs_data {

GridGeneratorINSIM::GridGeneratorINSIM(INSIMMeshConfig const & config, std::vector<WellConfig> const & wells)
    : _config(config), _wells(wells)
{
  if ( _wells.size() < 2 )
    throw std::invalid_argument("Cannot build INSIM grid with less than two wells");

  for (auto const & well : _wells) {
    if ( well.coordinates.size() == 1 )
      setup_simple_well_(well);
    else
      setup_complex_well_(well);
  }

  if (_config.padding_fraction <= 0.)
    throw std::invalid_argument("Invalid padding fraction for INSIM grid generator");

  generate_bounding_box_();
  extend_bounding_box_();
}

void GridGeneratorINSIM::extend_bounding_box_()
{
  // double const margin = _config.padding_fraction * find_characteristic_length_();
  // bbox_min -= margin;
  // bbox_max += margin;
  auto const c = _bbox->center();
  auto & vertices = _bbox->get_points();
  double const characteristic_length = c.distance( vertices[0] );
  // double const padding = _config.padding_fraction * characteristic_length;
  // We need to make the bounding box a perfect cube; otherwise, GMsh may generate really small cells
  angem::Point<3,double> dims{
    vertices[0].distance(vertices[1]),
    vertices[0].distance(vertices[3]),
    vertices[0].distance(vertices[4]),
  };
  // determine thickest dimension
  double maxdim = -1;
  for (size_t i = 0; i < 3; ++i)
    maxdim = std::max(maxdim, dims[i]);

  // maxdim += padding;

  // unify cube dimensions
  // move in x direction
  for (size_t v : {0, 3, 7, 4})
    vertices[v][0] -= 0.5*(maxdim - dims[0]);
  for (size_t v : {1, 5, 6, 2})
    vertices[v][0] += 0.5*(maxdim - dims[0]);
  // move in y direction
  for (size_t v : {0, 1, 5, 4})
    vertices[v][1] -= 0.5*(maxdim - dims[1]);
  for (size_t v : {2, 3, 7, 6})
    vertices[v][1] += 0.5*(maxdim - dims[1]);
  // move in z direction
  for (size_t v : {0, 1, 2, 3})
    vertices[v][2] -= 0.5*(maxdim - dims[2]);
  for (size_t v : {4, 5, 6, 7})
    vertices[v][2] += 0.5*(maxdim - dims[2]);
}

void GridGeneratorINSIM::assign_cell_labels_(mesh::Mesh & grid) const
{
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    cell->set_marker(_config.cell_label);
}


void GridGeneratorINSIM::generate_bounding_box_()
{
  double const upper = std::numeric_limits<double>::max();
  double const lower = std::numeric_limits<double>::lowest();
  angem::Point<3,double> bbox_min = {upper, upper, upper};
  angem::Point<3,double> bbox_max = {lower, lower, lower};
  for (auto const & v : _vertices)
    for (size_t i = 0; i < 3; ++i) {
      bbox_min[i] = std::min( bbox_min[i], v[i] );
      bbox_max[i] = std::max( bbox_max[i], v[i] );
    }

  double const padding = _config.padding_fraction * bbox_min.distance( bbox_max );
  bbox_min -= padding;
  bbox_max += padding;

  auto const delta = bbox_max - bbox_min;

  std::vector<angem::Point<3,double>> verts(8);  // hex has 8 vertices
  std::fill( verts.begin(), verts.begin() + 4, bbox_min );

  // look at the schematics in angem/Hexahedron.hpp to understand the order
  // to understand vertex numbering
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
  angem::Point<3,double> const c = angem::compute_center_mass( _vertices );

  // compute average distance between c and vertices
  double dist = 0.f;
  for (auto const & v : _vertices)
    dist += v.distance(c);

  return dist / _vertices.size();
}


GridGeneratorINSIM::operator mesh::Mesh() const
{
  mesh::Mesh grid;

#ifdef WITH_GMSH
  GmshInterface::initialize_gmsh(/*verbose = */ false);
  GmshInterface::build_triangulation_embedded_points(*_bbox, _vertices, grid);
  GmshInterface::finalize_gmsh();
#endif
  // mesh::IO::VTKWriter::write_geometry(grid, "test.vtk");

  assign_cell_labels_(grid);

  return grid;
}

void GridGeneratorINSIM::setup_simple_well_(WellConfig const & conf)
{
  _vertices.push_back(conf.coordinates[0]);
}

void GridGeneratorINSIM::setup_complex_well_(WellConfig const & conf)
{
  assert( conf.perforated.size() == conf.coordinates.size() - 1 );

  for (size_t i = 0; i < conf.perforated.size(); ++i)
    if ( conf.perforated[i] ) {
      auto const vertex = 0.5 * (conf.coordinates[i] + conf.coordinates[i+1]);
      _vertices.push_back(vertex);
    }
}

}  // end namespace gprs_data
