#include "GmshInterface.hpp"
#include <cstdlib> // atoi
#include <numeric>  // provides std;:iota

namespace gprs_data
{

using angem::VTK_ID;
using Point = angem::Point<3,double>;

// map gmsh element id to vtk id
std::vector<int> msh_id_to_vtk_id = {
    VTK_ID::InvalidElementID,     // [0]  does not exist
    VTK_ID::LineID,               // [1]  2-node edge
    VTK_ID::TriangleID,           // [2]  3-node triangle
    VTK_ID::QuadrangleID,         // [3]  4-node quadrangle
    VTK_ID::TetrahedronID,        // [4]  4-node tetrahedron
    VTK_ID::HexahedronID,         // [5]  8-node hexahedron
    VTK_ID::WedgeID,              // [6]  6-node prism
    VTK_ID::PyramidID,            // [7]  5-node pyramid
    VTK_ID::QuadraticEdgeID,      // [8]  3-node second order line
    VTK_ID::QuadraticTriangleID,  // [9]  6-node second order triangle
    VTK_ID::QuadraticQuadID,      // [10] 9-node second orderer quadrangle
    VTK_ID::QuadraticTetraID,     // [11] 10-node second order tetrahedron
    VTK_ID::InvalidElementID,     // [12] 27-node second order hexahedron
    VTK_ID::InvalidElementID,     // [13] 18-node second order prism
    VTK_ID::InvalidElementID,     // [14] 14-node second order pyramid
    VTK_ID::VertexID,             // [15] 1-node point
    VTK_ID::InvalidElementID,     // [16] 8-node second order quadrangle
    VTK_ID::QudraticHexahedronID, // [17] 20-node second order hexahedron
    VTK_ID::QuadraticWedgeID,     // [18] 15-node second order prism
    VTK_ID::QuadraticPyramidID,   // [19] 13-node second order pyramid
    VTK_ID::InvalidElementID,     // [20] 9-node third order incomplete triangle
    VTK_ID::InvalidElementID,     // [21] 10-node third order triangle
    VTK_ID::InvalidElementID,     // [22] 12-node fourth order incomplete triangle
    VTK_ID::InvalidElementID,     // [23] 15-node fourth order triangle
    VTK_ID::InvalidElementID,     // [24] 15-node fifth order incomplete triangle
    VTK_ID::InvalidElementID,     // [25] 21-node fifth order complete triangle
    VTK_ID::InvalidElementID,     // [26] 4-node third order edge
    VTK_ID::InvalidElementID,     // [27] 5-node fourth order edge
    VTK_ID::InvalidElementID,     // [28] 6-node fifth order edge
    VTK_ID::InvalidElementID,     // [29] 20-node third order tetrahedron
    VTK_ID::InvalidElementID,     // [30] 35-node fourth order tetrahedron
    VTK_ID::InvalidElementID      // [31] 56-node fifth order tetrahedron
};

enum GmshElementType
{
  node,
  edge,
  face,
  cell,
  invalid_element
};

std::vector<GmshElementType> gmsh_element_types = {
 GmshElementType::invalid_element, // [0]  does not exist
 GmshElementType::edge,            // [1]  2-node edge
 GmshElementType::face,            // [2]  3-node triangle
 GmshElementType::face,            // [3]  4-node quadrangle
 GmshElementType::cell,            // [4]  4-node tetrahedron
 GmshElementType::cell,            // [5]  8-node hexahedron
 GmshElementType::cell,            // [6]  6-node prism
 GmshElementType::cell,            // [7]  5-node pyramid
 GmshElementType::invalid_element, // [8]  3-node second order line
 GmshElementType::invalid_element, // [9]  6-node second order triangle
 GmshElementType::invalid_element, // [10] 9-node second orderer quadrangle
 GmshElementType::invalid_element, // [11] 10-node second order tetrahedron
 GmshElementType::invalid_element, // [12] 27-node second order hexahedron
 GmshElementType::invalid_element, // [13] 18-node second order prism
 GmshElementType::invalid_element, // [14] 14-node second order pyramid
 GmshElementType::node,            // [15] 1-node point
 GmshElementType::invalid_element, // [16] 8-node second order quadrangle
 GmshElementType::invalid_element, // [17] 20-node second order hexahedron
 GmshElementType::invalid_element, // [18] 15-node second order prism
 GmshElementType::invalid_element, // [19] 13-node second order pyramid
 GmshElementType::invalid_element, // [20] 9-node third order incomplete triangle
 GmshElementType::invalid_element, // [21] 10-node third order triangle
 GmshElementType::invalid_element, // [22] 12-node fourth order incomplete triangle
 GmshElementType::invalid_element, // [23] 15-node fourth order triangle
 GmshElementType::invalid_element, // [24] 15-node fifth order incomplete triangle
 GmshElementType::invalid_element, // [25] 21-node fifth order complete triangle
 GmshElementType::invalid_element, // [26] 4-node third order edge
 GmshElementType::invalid_element, // [27] 5-node fourth order edge
 GmshElementType::invalid_element, // [28] 6-node fifth order edge
 GmshElementType::invalid_element, // [29] 20-node third order tetrahedron
 GmshElementType::invalid_element, // [30] 35-node fourth order tetrahedron
 GmshElementType::invalid_element  // [31] 56-node fifth order tetrahedron
};

std::vector<size_t> gmsh_element_nvertices = {
 0,  // [0]  does not exist
 2,  // [1]  2-node edge
 3,  // [2]  3-node triangle
 4,  // [3]  4-node quadrangle
 4,  // [4]  4-node tetrahedron
 8,  // [5]  8-node hexahedron
 6,  // [6]  6-node prism
 5,  // [7]  5-node pyramid
 3,  // [8]  3-node second order line
 6,  // [9]  6-node second order triangle
 9,  // [10] 9-node second orderer quadrangle
 10, // [11] 10-node second order tetrahedron
 27, // [12] 27-node second order hexahedron
 18, // [13] 18-node second order prism
 14, // [14] 14-node second order pyramid
 1,  // [15] 1-node point
 8,  // [16] 8-node second order quadrangle
 20, // [17] 20-node second order hexahedron
 15, // [18] 15-node second order prism
 13, // [19] 13-node second order pyramid
 20, // [20] 9-node third order incomplete triangle
 10, // [21] 10-node third order triangle
 22, // [22] 12-node fourth order incomplete triangle
 15, // [23] 15-node fourth order triangle
 15, // [24] 15-node fifth order incomplete triangle
 21, // [25] 21-node fifth order complete triangle
 4,  // [26] 4-node third order edge
 5,  // [27] 5-node fourth order edge
 6,  // [28] 6-node fifth order edge
 20, // [29] 20-node third order tetrahedron
 35, // [30] 35-node fourth order tetrahedron
 56  // [31] 56-node fifth order tetrahedron
};


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

  // std::size_t n_cells, n_faces;
  std::string line;
  std::vector<std::string> strings;

  // const std::vector<int> polyhedras = {4, 5, 6, 11, 17, 18};
  // const std::vector<int> polygons = {2, 3, 8, 16};
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
    if ((element_type < 0) || (element_type > msh_id_to_vtk_id.size()))
      throw std::invalid_argument("Wrong vtk id type");
    const int vtk_id = msh_id_to_vtk_id[element_type];

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
        mesh.insert_face(ivertices, vtk_id, marker);
        break;
      case GmshElementType::cell:
        mesh.insert_cell(ivertices, vtk_id, marker);
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

    if ((element_type < 0) || (element_type > msh_id_to_vtk_id.size()))
      throw std::invalid_argument("Wrong vtk id type");
    const int vtk_id = msh_id_to_vtk_id[element_type];

    const int n_element_vertices = gmsh_element_nvertices[element_type];
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
        mesh.insert_face(vertices, vtk_id, physical_tag);
      if (entity_dim == 3)  // cells
        mesh.insert_cell(vertices, vtk_id, physical_tag);

      element++;
    }
  }

  std::cout << "n_cells = " << mesh.n_cells() << std::endl;
}

#ifdef WITH_GMSH

void GmshInterface::build_triangulation(const mesh::Cell & cell)
{
  const auto poly = cell.polyhedron();
  build_triangulation_(*poly);
}

void GmshInterface::build_triangulation_(const angem::Polyhedron<double> & cell)
{
  gmsh::option::setNumber("General.Terminal", 1);

  const double discr_element_size = 0.1 * compute_element_size_(cell);
  std::cout << "discr_element_size = " << discr_element_size << std::endl;
  // build points
  const std::vector<Point> & vertices = cell.get_points();
  std::cout << "n_vertices = " << vertices.size() << std::endl;
  for (std::size_t i=0; i < vertices.size(); ++i)
  {
    const Point & vertex = vertices[i];
    gmsh::model::geo::addPoint(vertex.x(), vertex.y(), vertex.z(),
                               discr_element_size, /*tag = */ i);
  }

  // build lines (edges)
  const auto edges = cell.get_edges();
  for (std::size_t i=0; i<edges.size(); ++i)
  {
    const std::pair<size_t,size_t> & edge = edges[i];
    gmsh::model::geo::addLine(edge.first, edge.second, i);
    // gmsh::model::addPhysicalGroup(1, {i}, i);
  }

  // build faces
  const std::vector<std::vector<std::size_t>> & faces = cell.get_faces();
  std::vector<angem::Polygon<double>> face_polys = cell.get_face_polygons();
  for (std::size_t i=0; i<face_polys.size(); ++i)
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
      assert(it_edge != edges.end());;
      const int edge_marker = static_cast<int>(std::distance(edges.begin(), it_edge));
      // figure out the sign of the tage
      if (cell_edge_ordered.first == cell_edge_unordered.first)
        edge_markers.push_back(edge_marker);
      else  // inverse orientation
        edge_markers.push_back(-edge_marker);
    }
    // create line loop and surface
    // NOTE: curve and surface loop must start from 1, otherwise gmsh
    // throws an error, ergo i+1
    gmsh::model::geo::addCurveLoop(edge_markers, static_cast<int>(i+1));
    gmsh::model::geo::addPlaneSurface({static_cast<int>(i+1)}, static_cast<int>(i+1));
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

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);
  gmsh::write("cell.msh");
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

void GmshInterface::initialize_gmsh()
{
  gmsh::initialize();
}

void GmshInterface::finalize_gmsh()
{
  gmsh::finalize();
}

#else

void GmshInterface::build_grid(const mesh::Cell & cell)
{
  throw std::invalid_argument("Gmsh is not linked. Gridding options not available");
}

void GmshInterface::initialize_gmsh()
{}

void GmshInterface::finalize_gmsh()
{}

#endif
}  // end namespace gprs_data
