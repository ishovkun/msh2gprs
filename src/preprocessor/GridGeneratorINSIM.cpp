#include "GridGeneratorINSIM.hpp"
#include "Well.hpp"                                       // provides Well
#include "gmsh_interface/GmshInterface.hpp"
#ifdef WITH_GMSH
#include <gmsh.h>
#endif

namespace gprs_data {

GridGeneratorINSIM::GridGeneratorINSIM(INSIMMeshConfig const & config, std::vector<WellConfig> const & wells)
    : _config(config), _wells(wells)
{
  if ( _wells.size() < 2 )
    throw std::invalid_argument("Cannot build INSIM grid with less than two wells");

  for (auto const & well : _wells) {
    std::cout << "well.name = " << well.name << std::endl;
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
  for (auto const & v : _vertices)
    for (size_t i = 0; i < 3; ++i) {
      bbox_min[i] = std::min( bbox_min[i], v[i] );
      bbox_max[i] = std::max( bbox_max[i], v[i] );
    }

  bbox_min -= margin;
  bbox_max += margin;

  auto const delta = bbox_min - bbox_max;

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
  angem::Point<3,double> const c = angem::compute_center_mass( _vertices );

  // compute average distance between c and vertices
  double dist = 0.f;
  for (auto const & v : _vertices)
    dist += v.distance(c);

  return dist / _vertices.size();
}


GridGeneratorINSIM::operator mesh::Mesh() const
{
#ifdef WITH_GMSH
  GmshInterface::initialize_gmsh(/*verbose = */ true);
  GmshInterface::build_triangulation_embedded_points(*_bbox, _vertices);
  // gmsh::initialize();
  // gmsh::option::setNumber("General.Terminal", 0);  // 0 shuts up gmsh logging
  // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  // gmsh::model::add("cell1");
  // gmsh::option::setNumber("Mesh.SaveAll", 1);

  // // outer points
  // const std::vector<angem::Point<3,double>> & vertices = _bbox->get_points();
  // double const characteristic_size = vertices.front().distance(vertices.back());
  // for (size_t i = 0; i < vertices.size(); ++i) {
  //   gmsh::model::geo::addPoint( vertices[i].x(), vertices[i].y(), vertices[i].z(),
  //                               characteristic_size, /*tag = */ i+1 );
  //   gmsh::model::addPhysicalGroup(0, {static_cast<int>(i+1)}, i+1);
  // }

  // // outer lines (edges)
  // const auto edges = _bbox.get_edges();
  // for (size_t i=0; i<edges.size(); ++i)
  // {
  //   const std::pair<size_t,size_t> & edge = edges[i];
  //   gmsh::model::geo::addLine(edge.first+1, edge.second+1, i+1);
  //   gmsh::model::addPhysicalGroup(1, {static_cast<int>(i+1)}, i+1);
  // }

  // // add embedded points (well vertices)
  // // for (auto const & v : _vertices) {
  //   gmsh::model::mesh::embed(3, {1,2}, 2, 3);
  // // }
  //     // Embed the model entities of dimension `dim' and tags `tags' in the
  //     // (`inDim', `inTag') model entity. The dimension `dim' can 0, 1 or 2 and
  //     // must be strictly smaller than `inDim', which must be either 2 or 3. The
  //     // embedded entities should not intersect each other or be part of the
  //     // boundary of the entity `inTag', whose mesh will conform to the mesh of the
  //     // embedded entities. With the OpenCASCADE kernel, if the `fragment'
  //     // operation is applied to entities of different dimensions, the lower
  //     // dimensional entities will be automatically embedded in the higher
  //     // dimensional entities if they are not on their boundary.
  //     // GMSH_API void embed(const int dim,
  //     //                     const std::vector<int> & tags,
  //     //                     const int inDim,
  //     //                     const int inTag);

  // gmsh::model::geo::synchronize();
  // gmsh::model::mesh::generate(2);
  // gmsh::write("cell.msh");
  // gmsh::finalize();


#endif
  return mesh::Mesh();
}

void GridGeneratorINSIM::setup_simple_well_(WellConfig const & conf)
{
  assert( _config.minimum_thickness > 0 );
  auto const & v1 = conf.coordinates[0];
  _vertices.push_back( v1 );
  auto v2 = v1;
  v2[2] -= _config.minimum_thickness;  // shift down in z direction
  _vertices.push_back( v2 );
}

}  // end namespace gprs_data
