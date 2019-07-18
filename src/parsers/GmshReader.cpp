#include <GmshReader.hpp>
#include <angem/PolyhedronFactory.hpp>
#include <fstream>
#include <sstream>      // std::stringstream
#include <ios>
#include <iterator>
#include <cstdlib> // atoi


namespace Parsers
{

// typedef std::unordered_map<int,int> MapIntInt;
// static int get_vtk_index(const int gmsh_index) {return map_gmsh_vtk[gmsh_index];}
GmshReader::MapIntInt GmshReader::map_gmsh_vtk = {
  // polyhedras
  {4, 10},  // tetra4
  {5, 12},  // prism8
  {6, 13},  // prism6
  {11, 24}, // tetra10
  {17, 25}, // prism20
  {18, 26}, // prism15
  // polygons
  {2, 5}, // trgle3
  {3, 9}, // quad4
  {9, 22}, // trgle6
  {16, 23}, // quad8
};


GmshReader::MapIntInt GmshReader::map_vtk_element_size = {
  // polyhedras
  {10, 4},  // tetra4
  {12, 8},  // prism8
  {13, 6},  // prism6
  {24, 10}, // tetra10
  {25, 20}, // prism20
  {26, 15}, // prism15
  // polygons
  {5, 3}, // trgle3
  {9, 4}, // quad4
  {22, 6}, // trgle6
  {23, 8}, // quad8
};


void GmshReader::read_input(const std::string & filename,
                            mesh::Mesh        & mesh)
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
    read_gmsh2_input(mesh_file, mesh);
  else if (version > 4.0)  // works for 4.1
    read_gmsh4_input(mesh_file, mesh);
  else
  {
    mesh_file.close();
    throw std::out_of_range("cannot read this msh file "
                            "(version "+ std::to_string(version) +")");
  }

  mesh_file.close();
}


void GmshReader::read_gmsh2_input(std::fstream & mesh_file,
                                  mesh::Mesh   & mesh)
{
  std::string entry;

  // read until nodes
  while(entry != "$Nodes")
    mesh_file >> entry;

  std::size_t n_vertices;
  mesh_file >> n_vertices;

  std::cout << "\tn_vertices = " << n_vertices << std::endl;
  mesh.vertices.points.reserve(n_vertices);

  const int dim = 3;
  for (std::size_t ivertex=0; ivertex<n_vertices; ++ivertex)
  {
    double value;
    mesh_file >> value;  // skip vertex index
    angem::Point<3,double> vertex;
    for (int d=0; d<dim; ++d)
      mesh_file >> vertex[d];

    // warning if duplicates
    const std::size_t new_ind = mesh.vertices.insert(vertex);
    if (new_ind != mesh.vertices.points.size() - 1)
      std::cout << "WARNING: duplicate entry in gmsh file "
                << value << "\t"
                << vertex
                << std::endl;
    // mesh.vertices.points.push_back(vertex);
  }
  // endonodes

  // read until elements
  while(entry != "$Elements")
    mesh_file >> entry;

  std::size_t n_elements;
  mesh_file >> n_elements;

  std::cout << "\tn_elements = " << n_elements << std::endl;
  mesh.cells.reserve(n_elements);

  // std::size_t n_cells, n_faces;
  std::string line;
  std::vector<std::string> strings;

  const std::vector<int> polyhedras = {4, 5, 6, 11, 17, 18};
  const std::vector<int> polygons = {2, 3, 8, 16};
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
    const int vtk_id = get_vtk_index(element_type);
    const int marker = std::atoi(tokens[3].c_str());
    std::vector<std::size_t> ivertices(tokens.size() - vert_shift);
    for (int j=0; j<tokens.size() - vert_shift; ++j)
      ivertices[j] = std::atoi(tokens[j+vert_shift].c_str()) - 1;

    // 3D element
    if (std::find(begin(polyhedras), end(polyhedras), element_type) != polyhedras.end())
    {
      // const auto polyhedron = angem::PolyhedronFactory::create
      //     (mesh.vertices.points, ivertices, vtk_id);
      // mesh.insert(polyhedron, marker);
      mesh.insert_cell(ivertices, vtk_id, marker);
    }
    // 2D element
    else if (std::find(begin(polygons), end(polygons), element_type) !=
             polygons.end())
    {
      // const angem::Polygon<double> poly(mesh.vertices.points, ivertices);
      // mesh.insert(poly, marker);
      mesh.insert_face(ivertices, vtk_id, marker);
    }
    else
      throw std::out_of_range("unknown element type");
  }
}


void GmshReader::read_gmsh4_input(std::fstream & mesh_file,
                                  mesh::Mesh   & mesh)
{
  std::string entry;

  //  read until nodes
  while(entry != "$Nodes")
    mesh_file >> entry;

  // n_entities whatever that is
  mesh_file >> entry;

  // read number of vertices
  std::size_t n_vertices;
  mesh_file >> n_vertices;
  std::cout << "\tn_vertices = " << n_vertices << std::endl;
  mesh.vertices.points.reserve(n_vertices);

  const int dim = 3;
  for (std::size_t ivertex=0; ivertex<n_vertices; ++ivertex)
  {
    double value;
    mesh_file >> value;  // skip vertex index
    angem::Point<3,double> vertex;
    for (int d=0; d<dim; ++d)
      mesh_file >> vertex[d];

    // std::cout << "adding vertex " << vertex << std::endl;
    //  warning if duplicates
    const std::size_t new_ind = mesh.vertices.insert(vertex);
    if (new_ind != mesh.vertices.points.size() - 1)
      std::cout << "WARNING: duplicate entry in gmsh file "
                << value << "\t"
                << vertex
                << std::endl;
    // mesh.vertices.points.push_back(vertex);
  } // endonodes

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

  std::string line;
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
    const int element_type = std::atoi(tokens[2].c_str());
    const int vtk_id = get_vtk_index(element_type);
    const int n_element_vertices = map_vtk_element_size[vtk_id];

    std::size_t n_elements_in_block;
    std::stringstream sstream(tokens[3]);
    sstream >> n_elements_in_block;
    std::vector<std::size_t> vertices(n_element_vertices);

    if (entity_dim == 3)
      mesh.cells.reserve(mesh.cells.capacity() + n_elements);

    for (size_t i=0; i< n_elements_in_block; i++)
    {
      // line 3: element tag, element nodes
      getline(mesh_file, line);
      if (entity_dim == 2)  // fuck 2d elements
      {}
      if (entity_dim == 3)  // fuck 2d elements
      {
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};
        static const int vert_shift = 1;  // index of first vertex in line
        for (int j=vert_shift; j<tokens.size(); ++j)
          vertices[j] = std::atoi(tokens[j].c_str()) - 1;

        mesh.insert_cell(vertices, vtk_id, 0);
      }

      element++;
    }
  }

  std::cout << "n_cells = " << mesh.n_cells() << std::endl;

  throw std::invalid_argument("reading gmsh v.4 not implemented yet");
}

}
