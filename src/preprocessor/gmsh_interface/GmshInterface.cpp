#include "GmshInterface.hpp"
#include <cstdlib> // atoi
#include <numeric>  // provides std;:iota
#include <algorithm> // max_element
#include <sstream>   // isstream
#include <iterator>  // std::istream_iterator

namespace gprs_data
{

using Point = angem::Point<3,double>;

void GmshInterface::read_msh(const std::string & filename, mesh::Mesh & mesh)
{
  std::fstream mesh_file;
  mesh_file.open(filename.c_str(), std::fstream::in);
  if (!mesh_file)
    throw std::out_of_range(filename + " does not exist");

  std::string line;
  // read mesh format
  std::getline(mesh_file, line);
  float version;
  mesh_file >> version;
  std::cout << "Gmsh file version: " << version << std::endl;

  if (version >= 2 and version < 3)
    read_msh_v2_(mesh_file, mesh);
  else if (version > 4.0)  // works for 4.1
  {
    std::cout << "warning, gmsh 4 files not tested" << std::endl;
    read_msh_v4_(mesh_file, mesh);
  }
  else
  {
    mesh_file.close();
    throw std::out_of_range("cannot read this msh file "
                            "(version "+ std::to_string(version) +")");
  }

  mesh_file.close();
}

void GmshInterface::read_msh_v2_(std::fstream & mesh_file, mesh::Mesh & mesh)
{
  std::string entry;

  // read until nodes
  while(entry != "$Nodes")
    mesh_file >> entry;

  std::size_t n_vertices;
  mesh_file >> n_vertices;

  std::cout << "\tn_vertices = " << n_vertices << std::endl;
  auto & vertices= mesh.vertices();
  mesh.vertices().reserve(n_vertices);
  angem::PointSet<3,double> pset;

  const int dim = 3;
  for (std::size_t ivertex=0; ivertex<n_vertices; ++ivertex)
  {
    double value;
    mesh_file >> value;  // skip vertex index
    angem::Point<3,double> vertex;
    for (int d=0; d<dim; ++d)
      mesh_file >> vertex[d];
    size_t ins = pset.insert(vertex);
    if (ins != pset.size() - 1)
    {
      std::cout << "Error. duplicated vertex" << std::endl;
      exit(0);
    }

    vertices.push_back(vertex);
  }
  // endonodes

  // read until elements
  while(entry != "$Elements")
    mesh_file >> entry;

  std::size_t n_elements;
  mesh_file >> n_elements;

  std::cout << "\tn_elements = " << n_elements << std::endl;
  mesh.cells().reserve(n_elements);

  std::string line;
  std::vector<std::string> strings;

  for (std::size_t i=0; i<n_elements; ++i)
  {
    getline(mesh_file, line);
    std::istringstream iss(line);
    std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
    if (tokens.empty())
    {
      i--;
      continue;
    }

    // tokens:
    // first - element number
    // second - type
    // third - ?
    // fourth - marker
    // fifth - ?
    // later - vertex indices
    static const int vert_shift = 5;  // index of first vertex in line
    const int element_type = std::atoi(tokens[1].c_str());
    const int marker = std::atoi(tokens[3].c_str());
    std::vector<std::size_t> ivertices(tokens.size() - vert_shift);
    for (int j=0; j<tokens.size() - vert_shift; ++j)
      ivertices[j] = std::atoi(tokens[j+vert_shift].c_str()) - 1;

    // 3D element
    if (element_type >= gmsh_element_types.size() or element_type < 0)
    {
      const std::string msg = "Unknown element type: "  + std::to_string(element_type);
      throw std::out_of_range(msg);
    }
    switch (gmsh_element_types[element_type])
    {
      case GmshElementType::node:
        continue;
      case GmshElementType::edge:
        continue;
      case GmshElementType::face:
        mesh.insert_face(ivertices, get_vtk_id(element_type), marker);
        break;
      case GmshElementType::cell:
        mesh.insert_cell(ivertices, get_vtk_id(element_type), marker);
        break;
      default:
        const std::string msg = "Unknown element type: "  + std::to_string(element_type);
        throw std::out_of_range(msg);
    }
  }
}

void GmshInterface::read_msh_v4_(std::fstream & mesh_file, mesh::Mesh & mesh)
{
  std::string entry;

  //  read until entities
  while(entry != "$Entities") mesh_file >> entry;

  /* entities describe tags that we need in some cases
   * We primarily interested in surface tags (bc's, fracs)
   * and volume tags (reservoir subdomains)
   *
   * numPoints(size_t) numCurves(size_t)
   * numSurfaces(size_t) numVolumes(size_t)
   *
   * pointTag(int) X(double) Y(double) Z(double)
   *      numPhysicalTags(size_t) physicalTag(int) ...
   *
   * curveTag(int) minX(double) minY(double) minZ(double)
   *     maxX(double) maxY(double) maxZ(double)
   *     numPhysicalTags(size_t) physicalTag(int) ...
   *     numBoundingPoints(size_t) pointTag(int) ...  //
   *
   * surfaceTag(int) minX(double) minY(double) minZ(double)
   *     maxX(double) maxY(double) maxZ(double)
   *     numPhysicalTags(size_t) physicalTag(int) ...
   *     numBoundingCurves(size_t) curveTag(int) ...
   *
   * volumeTag(int) minX(double) minY(double) minZ(double)
   *     maxX(double) maxY(double) maxZ(double)
   *     numPhysicalTags(size_t) physicalTag(int) ...
   *     numBoundngSurfaces(size_t) surfaceTag(int) ...
   */

  std::string line;
  std::vector<std::string> tokens;
  while (tokens.empty())
  {
    getline(mesh_file, line);
    std::istringstream iss(line);
    tokens = {std::istream_iterator<std::string>{iss},
              std::istream_iterator<std::string>{}};
  }
  if (tokens.size() != 4)
  {
    std::cout << line << std::endl;
    throw std::invalid_argument("invalid mesh format");
  }

  std::size_t n_points, n_curves, n_surfaces, n_volumes;
  {std::stringstream st(tokens[0]); st >> n_points;}
  {std::stringstream st(tokens[1]); st >> n_curves;}
  {std::stringstream st(tokens[2]); st >> n_surfaces;}
  {std::stringstream st(tokens[3]); st >> n_volumes;}

  std::cout << "n_points = " << n_points << std::endl;
  std::cout << "n_curves = " << n_curves << std::endl;
  std::cout << "n_physical_surfaces = " << n_surfaces << std::endl;
  std::cout << "n_subdomains = " << n_volumes << std::endl;

  for (std::size_t i = 0; i < n_points; i++)  // skip points
  {
    getline(mesh_file, line);
    // std::cout << line << std::endl;
  }
  for (std::size_t i = 0; i < n_curves; i++)  // skip curves
  {
    getline(mesh_file, line);
    // std::cout << line << std::endl;
  }

  //  read surfaces
  std::unordered_map<int, int> surface_tags;
  for (std::size_t i = 0; i < n_surfaces; i++)
  {
   /* surfaceTag(int) minX(double) minY(double) minZ(double)
   *     maxX(double) maxY(double) maxZ(double)
   *     numPhysicalTags(size_t) physicalTag(int) ...
   *     numBoundingCurves(size_t) curveTag(int) ...
   */
    getline(mesh_file, line);
    std::istringstream iss(line);
    std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
    const int entity = std::atoi(tokens[0].c_str());
    const int n_physical_tags = std::atoi(tokens[7].c_str());
    const int tag = (n_physical_tags == 1) ?
        std::atoi(tokens[8].c_str()) : mesh::constants::default_face_marker;
    if (n_physical_tags > 1)
    {
      std::cout << "error in line: " << line << std::endl;
      std::cout << "n_physical_tags = " << n_physical_tags << std::endl;
      throw std::invalid_argument("more than one physical tag per surface not supported");
    }

    surface_tags.insert({entity, tag});
  }

  // read volumes
  std::unordered_map<int, int> volume_tags;
  for (std::size_t i = 0; i < n_volumes; i++)
  {
    getline(mesh_file, line);
    std::istringstream iss(line);
    std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
    const int entity = std::atoi(tokens[0].c_str());
    const int n_physical_tags = std::atoi(tokens[7].c_str());
    if (n_physical_tags != 1)
      throw std::invalid_argument("more than one physical tag per volume not supported");
    const int tag = std::atoi(tokens[8].c_str());
    volume_tags.insert({entity, tag});
  }

  //  read until nodes
  while(entry != "$Nodes")
    mesh_file >> entry;

  /* Vertex data as follows:
   * numEntityBlocks(size_t) numNodes(size_t)
   * minNodeTag(size_t) maxNodeTag(size_t)
   * entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
   * nodeTag(size_t)
   * ...
   * x(double) y(double) z(double) */

  mesh_file >> entry;  // numEntityBlocks
  std::size_t n_vertices; mesh_file >> n_vertices;  // numNodes
  mesh_file >> entry;  // minNodeTag
  mesh_file >> entry;  // maxNodeTag

  std::cout << "\tn_vertices = " << n_vertices << std::endl;
  mesh.vertices().reserve(n_vertices);

  size_t vertex = 0;
  while (vertex < n_vertices)
  {
    mesh_file >> entry;  // entityDim (0 for vertices)
    mesh_file >> entry;  // entityTag
    mesh_file >> entry;  // parametric
    std::size_t n_nodes_in_block;
    mesh_file >> n_nodes_in_block;  // numNodesInBlock(size_t)

    // skip node tags
    for (std::size_t j=0; j<n_nodes_in_block; ++j)
    {
      size_t node_tag;
      mesh_file >> node_tag;  // nodeTag
    }

    // read vertices
    for (std::size_t j=0; j<n_nodes_in_block; ++j)
    {
      static const int dim = 3;
      angem::Point<dim,double> coords;
      for (int d=0; d<dim; ++d) mesh_file >> coords[d];

      mesh.vertices().push_back(coords);
      vertex++;
    } // endonodes
  }


  //  read until elements
  while(entry != "$Elements")
    mesh_file >> entry;

  // line 1:
  //  numEntityBlocks(size_t) numElements(size_t)
  // minElementTag(size_t) maxElementTag(size_t)

  // number of entity blocks
  mesh_file >> entry;

  // number of elements
  std::size_t n_elements;
  mesh_file >> n_elements;

  std::cout << "\tn_elements = " << n_elements << std::endl;

  //  min element tag
  mesh_file >> entry;
  //  max element tag
  mesh_file >> entry;

  size_t element = 0;
  while (element < n_elements)
  {
    // line 2: block of data
    // entityDim(int) entityTag(int)
    // elementType(int; see below) numElementsInBlock(size_t)
    getline(mesh_file, line);
    std::istringstream iss(line);
    std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
    if (tokens.empty())
      continue;

    const int entity_dim = std::atoi(tokens[0].c_str());
    const int entity_tag = std::atoi(tokens[1].c_str());
    const int element_type = std::atoi(tokens[2].c_str());
    const int n_element_vertices = get_n_vertices(element_type);

    int physical_tag;
    if (entity_dim == 2)
      physical_tag = surface_tags[entity_tag];
    if (entity_dim == 3)
      physical_tag = volume_tags[entity_tag];

    std::size_t n_elements_in_block;
    std::stringstream sstream(tokens[3]);
    sstream >> n_elements_in_block;
    std::vector<std::size_t> vertices(n_element_vertices);

    if (entity_dim == 3)
      mesh.cells().reserve(mesh.cells().capacity() + n_elements);

    for (size_t i = 0; i < n_elements_in_block; i++)
    {
      // line 3: element tag, element nodes
      getline(mesh_file, line);
      std::istringstream iss(line);
      std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                      std::istream_iterator<std::string>{}};
      static const int vert_shift = 1;  // index of first vertex in line

      // fill node indices
      for (int j = vert_shift; j<tokens.size(); j++)
        vertices[j-vert_shift] = std::atoi(tokens[j].c_str()) - 1;

      if (entity_dim == 2)  // faces
        mesh.insert_face(vertices, get_vtk_id(element_type), physical_tag);
      if (entity_dim == 3)  // cells
        mesh.insert_cell(vertices, get_vtk_id(element_type), physical_tag);

      element++;
    }
  }

  std::cout << "n_cells = " << mesh.n_active_cells() << std::endl;
}

size_t GmshInterface::get_n_vertices(const int element_type)
{
  if ((element_type < 0) || (element_type > gmsh_element_nvertices.size()))
    throw std::invalid_argument("Wrong vtk id type");
  return gmsh_element_nvertices[element_type];
}

int GmshInterface::get_gmsh_element_id(const angem::VTK_ID vtk_id)
{
  const auto it = std::find( msh_id_to_vtk_id.begin(), msh_id_to_vtk_id.end(), vtk_id );
  if (it != msh_id_to_vtk_id.end())
    return std::distance(msh_id_to_vtk_id.begin(), it);
  else throw std::invalid_argument( "vtk id is not mapped" );
}


int GmshInterface::get_vtk_id(const int element_type)
{
  if ((element_type < 0) || (element_type > msh_id_to_vtk_id.size()))
    throw std::invalid_argument("Wrong vtk id type");
  return msh_id_to_vtk_id[element_type];
}

#ifdef WITH_GMSH

void GmshInterface::build_triangulation(const mesh::Cell & cell, const double n_vertices_on_edge)
{
  // const auto poly = cell.polyhedron();
  // build_triangulation_(*poly);
  std::vector<angem::Point<3,double>> coord = cell.vertex_coordinates();
  const auto verts = cell.vertices();
  std::vector<std::vector<size_t>> faces;
  for (auto face : cell.faces())
  {
    std::vector<size_t> poly_face;
    for (auto v : face->vertices())
    {
      const size_t idx = std::distance(verts.begin(), std::find( verts.begin(), verts.end(), v ));
      poly_face.push_back(idx);
    }
    faces.push_back(std::move(poly_face));
  }
  const angem::Polyhedron<double> poly(coord, faces);
  build_triangulation_(poly, n_vertices_on_edge);
}

std::vector<double>  compute_vertex_element_sizes(const std::vector<std::pair<size_t,size_t>> &edges,
                                                  const std::vector<angem::Point<3,double>> & coord)
{
  std::vector<double> result(coord.size(), std::numeric_limits<double>::max());
  for (const auto & edge : edges)
  {
    const double dist = coord[edge.first].distance(coord[ edge.second ]);
    result[edge.first] = std::min(result[ edge.first ], dist);
    result[edge.second] = std::min(result[ edge.second ], dist);
  }
  return result;
}

void GmshInterface::build_triangulation_(const angem::Polyhedron<double> & cell,
                                         const double n_vertices_on_edge)
{
  // gmsh::option::setNumber("General.Terminal", 1);
  gmsh::option::setNumber("General.Terminal", 0);  // 0 shuts up gmsh logging
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::model::add("cell1");

  // Remember that by default, if physical groups are defined, Gmsh will export
  // in the output mesh file only those elements that belong to at least one
  // physical group. To force Gmsh to save all elements, you can use
  //
  gmsh::option::setNumber("Mesh.SaveAll", 1);

  // const double discr_element_size = 0.4 * compute_element_size_(cell);
  // std::cout << "discr_element_size = " << discr_element_size << std::endl;
  // build points
  const std::vector<Point> & vertices = cell.get_points();
  const auto edges = cell.get_edges();
  const std::vector<double> element_sizes = compute_vertex_element_sizes(edges,vertices);

  for (size_t i=0; i < vertices.size(); ++i)
  {
    const Point & vertex = vertices[i];
    gmsh::model::geo::addPoint(vertex.x(), vertex.y(), vertex.z(),
                               element_sizes[i] / n_vertices_on_edge,
                               /*tag = */ i+1);
    gmsh::model::addPhysicalGroup(0, {static_cast<int>(i+1)}, i+1);
  }

  // build lines (edges)
  for (size_t i=0; i<edges.size(); ++i)
  {
    const std::pair<size_t,size_t> & edge = edges[i];
    gmsh::model::geo::addLine(edge.first+1, edge.second+1, i+1);
    gmsh::model::addPhysicalGroup(1, {static_cast<int>(i+1)}, i+1);
  }

  gmsh::model::geo::synchronize();
  for (size_t i=0; i<edges.size(); ++i)
    gmsh::model::mesh::setTransfiniteCurve(i+1, 5);

  // build faces
  const std::vector<std::vector<std::size_t>> & faces = cell.get_faces();
  std::vector<angem::Polygon<double>> face_polys = cell.get_face_polygons();
  for (size_t i=0; i<faces.size(); ++i)
  {
    const auto & face = faces[i];
    const auto & poly = face_polys[i];
    std::vector<int> edge_markers;
    for (const angem::Edge & face_edge : poly.get_edges())
    {
      const std::pair<size_t,size_t> cell_edge_ordered = {face[face_edge.first],
                                                         face[face_edge.second]};
      const std::pair<size_t,size_t> cell_edge_unordered = std::minmax(face[face_edge.first],
                                                                       face[face_edge.second]);
      // get edge marker
      const auto it_edge = std::find_if( edges.begin(), edges.end(),
                    [&cell_edge_unordered](const auto & it)->bool
                    {
                      return it.first == cell_edge_unordered.first &&
                             it.second == cell_edge_unordered.second;
                    });
      assert(it_edge != edges.end());
      const int edge_marker = static_cast<int>(std::distance(edges.begin(), it_edge));
      // figure out the sign of the edge
      if (cell_edge_ordered.first == cell_edge_unordered.first)
      {
        assert( cell_edge_ordered.second == cell_edge_unordered.second );
        edge_markers.push_back(edge_marker+1);
      }
      else  // inverse orientation
        edge_markers.push_back(-(edge_marker+1));
    }
    // create line loop and surface
    // NOTE: curve and surface loop must start from 1, otherwise gmsh
    // throws an error, ergo i+1
    // std::cout << "add Line loop " << i+1 << ":";
    // for (auto edge : edge_markers)
    //   std::cout << edge <<  " ";
    // std::cout << std::endl;
    gmsh::model::geo::addCurveLoop(edge_markers, static_cast<int>(i+1));
    // std::cout << "add plane surface " << i+1 << std::endl;
    gmsh::model::geo::addPlaneSurface({static_cast<int>(i+1)}, static_cast<int>(i+1));
    gmsh::model::addPhysicalGroup(2, {static_cast<int>(i+1)}, i+1);
  }

  // gmsh::model::geo::synchronize();
  // gmsh::model::mesh::generate(2);
  // gmsh::write("cell.msh");
  // gmsh::finalize();

  // create volume from surfaces
  std::vector<int> surfaces(faces.size());
  std::iota(surfaces.begin(), surfaces.end(), 1);
  gmsh::model::geo::addSurfaceLoop(surfaces, 1);
  gmsh::model::geo::addVolume({1}, 1);
  gmsh::model::addPhysicalGroup(3, {static_cast<int>(1)}, 1);

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);
  // gmsh::write("cell.msh");
  // // gmsh::finalize();
}

double GmshInterface::compute_element_size_(const angem::Polyhedron<double> & cell)
{
  double min_edge_length = std::numeric_limits<double>::max();
  const std::vector<Point> & vertices = cell.get_points();
  for (const auto & edge : cell.get_edges())
  {
    const double h = vertices[edge.first].distance(vertices[edge.second]);
    min_edge_length = std::min( min_edge_length, h);
  }
  return min_edge_length;
}

void GmshInterface::save_gmsh_grid(const std::string fname)
{
  gmsh::write(fname);
}


void GmshInterface::initialize_gmsh()
{
  gmsh::initialize();
}

void GmshInterface::finalize_gmsh()
{
  gmsh::finalize();
}

void GmshInterface::get_elements(std::vector<int> & element_types,
                                 std::vector<std::vector<std::size_t> > & element_tags,
                                 std::vector<std::vector<std::size_t> > & node_tags,
                                 const int dim,
                                 const int tag)
{
  gmsh::model::mesh::getElements(element_types, element_tags, node_tags, dim, tag);
}

void GmshInterface::build_triangulation(const mesh::Cell & cell, mesh::Mesh & grid,
                                        const double n_vertices_on_edge)
{
  build_triangulation(cell, n_vertices_on_edge);

  if (!grid.empty())
    throw std::invalid_argument("refuse to add gmsh elements to a non-empty grid");

  std::vector<size_t> node_tags;
  std::vector<double> node_coord, parametric_coord;
  gmsh::model::mesh::getNodes(node_tags, node_coord, parametric_coord, /*dim = */ 3,
                              /* tag */ -1, /*includeBoundary =*/ true,
                              /* return_parametric =  */ false);

  auto it = std::max_element(node_tags.begin(), node_tags.end());
  assert( it != node_tags.end() );
  const size_t max_node_tag = *it + 1;
  // find maximum node tag
  std::vector<size_t> vertex_numbering(max_node_tag, std::numeric_limits<size_t>::max());
  size_t iv = 0;
  for (const size_t node : node_tags)
    vertex_numbering[node]  = iv++;

  grid.vertices().reserve(node_tags.size());
  for (std::size_t i=0; i<node_tags.size(); ++i)
  {
    assert( node_coord[3*i+2] < node_coord.size() );
    angem::Point<3,double> vertex = { node_coord[3*i], node_coord[3*i+1], node_coord[3*i+2] };
    grid.vertices().push_back(vertex);
  }

  // insert cells
  insert_elements_(3, -1, vertex_numbering, grid);
  // insert face labels
  for (size_t parent_face=0; parent_face<cell.faces().size(); ++parent_face)
  {
    insert_elements_(2, parent_face+1, vertex_numbering, grid);
  }

}

void GmshInterface::insert_elements_(const int dim, const int tag,
                                     const std::vector<size_t> & vertex_numbering,
                                     mesh::Mesh & grid)
{
  std::vector<int> element_types;
  std::vector<std::vector<size_t> > element_tags;
  std::vector<std::vector<size_t> > element_node_tags;
  get_elements(element_types, element_tags, element_node_tags, dim, tag);
  for (size_t itype=0; itype<element_types.size(); ++itype)
  {
    const int element_type = element_types[itype];
    const int vtk_id = get_vtk_id(element_type);
    const size_t nv = get_n_vertices(element_type);
    // std::cout << "element_tags[itype].size() = " << element_tags[itype].size() << std::endl;

    for (size_t ie=0; ie < element_tags[itype].size(); ++ie)
    {
      // if (dim == 3) std::cout << "element_tags[itype] = " << element_tags[itype][ie] << std::endl;
      // get vertex indices from gmsh vertex tags
      std::vector<size_t> verts;
      verts.reserve(nv);
      for (size_t iv=0; iv<nv; ++iv)
      {
        const size_t vertex_tag = element_node_tags[itype][nv * ie + iv];
        const size_t vertex = vertex_numbering[vertex_tag];
        assert (vertex < grid.n_vertices() );
        verts.push_back(vertex);
      }

      // if (gmsh_element_types[element_type] == GmshElementType::face)
      //   std::cout << "face tag = " << element_tags[itype][ie] << std::endl;
      // if (gmsh_element_types[element_type] == GmshElementType::cell)
      //   std::cout << "cell tag = " << element_tags[itype][ie] << std::endl;

      // const int marker = (tag == -1) ? -1 : tag;
      switch (gmsh_element_types[element_type])
      {
        case GmshElementType::node:
          continue;
        case GmshElementType::edge:
          continue;
        case GmshElementType::face:
          grid.insert_face(verts, get_vtk_id(element_type), tag);
          break;
        case GmshElementType::cell:
          grid.insert_cell(verts, get_vtk_id(element_type), tag);
          break;
        default:
          const std::string msg = "Unknown element type: " + std::to_string(element_type);
          throw std::out_of_range(msg);
      }
    }
  }
}


#else

void GmshInterface::build_triangulation(const mesh::Cell & cell, const double n_vertices_on_edge)
{
  throw std::invalid_argument("Gmsh is not linked. Gridding options not available");
}

void GmshInterface::initialize_gmsh()
{}

void GmshInterface::finalize_gmsh()
{}

void GmshInterface::get_elements(std::vector<int> & elemen_types,
                                 std::vector<std::vector<std::size_t> > & element_tags,
                                 std::vector<std::vector<std::size_t> > & node_tags,
                                 const int dim,
                                 const int tag)
{}

static void get_mesh_from_model(mesh:: Mesh & grid)
{
  throw std::runtime_error("Gmsh has not been linked");
}

void GmshInterface::build_triangulation(const mesh::Cell & cell, mesh::Mesh & grid,
                                        const double n_vertices_on_edge)
{}


#endif
}  // end namespace gprs_data
