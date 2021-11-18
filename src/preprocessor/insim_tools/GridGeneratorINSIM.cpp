#include "GridGeneratorINSIM.hpp"
#include "Well.hpp"                                       // provides Well
#include "gmsh_interface/GmshInterface.hpp"
#include "MitchellBestCandidate.hpp"
#include "logger/Logger.hpp"
#include "mesh/io/VTKWriter.hpp"

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

  if (_config.padding_fraction < 0.)
    throw std::invalid_argument("Invalid padding fraction for INSIM grid generator");

  auto box = generate_bounding_box_();
  add_imaginary_wells_(box);
  // mesh::IO::VTKWriter::write(box, "box_orig.vtk");
  box = pad_bounding_box_(box);
  // mesh::IO::VTKWriter::write(box, "box_padded.vtk");
  box = extend_bounding_box_(box);
  // mesh::IO::VTKWriter::write(box, "box_extended.vtk");
  _bounding_box = std::make_unique<angem::Hexahedron<double>>(box);
}

void GridGeneratorINSIM::add_imaginary_wells_(angem::Hexahedron<double> const & box)
{
  if ( _config.n_imaginary_wells > 0 ) {
    logging::log() << "generating " << _config.n_imaginary_wells << " imaginary wells" << std::endl;
    MitchellBestCandidate mbc(_config);
    _vertices = mbc.generate_points(_vertices, box);
  }
}

angem::Hexahedron<double> GridGeneratorINSIM::extend_bounding_box_(angem::Hexahedron<double> const & original) const
{
  auto ans = original;
  auto const c = ans.center();
  auto & vertices = ans.get_points();
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

  return ans;
}

void GridGeneratorINSIM::assign_cell_labels_(mesh::Mesh & grid) const
{
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    cell->set_marker(_config.cell_label);
}


angem::Hexahedron<double> GridGeneratorINSIM::generate_bounding_box_() const
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
  // _bbox = std::make_unique<angem::Hexahedron<double>>(verts, indices);
  return angem::Hexahedron<double>(verts, indices);
}

GridGeneratorINSIM::operator mesh::Mesh() const
{
  mesh::Mesh grid;

#ifdef WITH_GMSH
  GmshInterface::initialize_gmsh(/*verbose = */ false);
  GmshInterface::build_triangulation_embedded_points( *_bounding_box, _vertices, grid );

  std::vector<size_t> node_tags;
  std::vector<double> node_coord, parametric_coord;
  // gmsh::model::mesh::getNodes(node_tags, node_coord, parametric_coord, -1,
  //                             /* tag */ -1, /* includeBoundary = */ true,
  //                             /* return_parametric =  */ false);
  // std::cout << "beer" << std::endl;
  // for (size_t i = 0; i < node_tags.size(); ++i)
  //   std::cout << node_tags[i]
  //             << " " << node_coord[3*i+0]
  //             << " " << node_coord[3*i+1]
  //             << " " << node_coord[3*i+2]
  //             << std::endl;

  GmshInterface::finalize_gmsh();
#endif
  // mesh::IO::VTKWriter::write_geometry(grid, "test.vtk");


  // exit(0);

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

angem::Hexahedron<double> GridGeneratorINSIM::pad_bounding_box_(angem::Hexahedron<double> const & original) const
{
  double const shift = original.radius() * _config.padding_fraction;

  auto ans = original;
  auto & vertices = ans.get_points();
  // move in x direction
  for (size_t v : {0, 3, 7, 4})
    vertices[v][0] -= shift;
  for (size_t v : {1, 5, 6, 2})
    vertices[v][0] += shift;
  // move in y direction
  for (size_t v : {0, 1, 5, 4})
    vertices[v][1] -= shift;
  for (size_t v : {2, 3, 7, 6})
    vertices[v][1] += shift;
  // move in z direction
  for (size_t v : {0, 1, 2, 3})
    vertices[v][2] -= shift;
  for (size_t v : {4, 5, 6, 7})
    vertices[v][2] += shift;

  return ans;
}

}  // end namespace gprs_data
