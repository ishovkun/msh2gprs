#pragma once

#include <mesh/Mesh.hpp>
#include <fstream>  // fstream

#ifdef WITH_GMSH
#include <gmsh.h>
#endif

namespace gprs_data
{

class GmshInterface
{
 public:
  // read .msh input file and save the result into grid
  static void read_msh(const std::string & filename, mesh::Mesh & grid);
  // run before building grid
  static void initialize_gmsh();
  // run after wrapping up the work with the grid
  static void finalize_gmsh();
  // build gmsh grid bounded by the faces of the cell
  static void build_triangulation(const mesh::Cell & cell);
  // Get the elements classified on the entity of dimension `dim' and tag
  // `tag'. If `tag' < 0, get the elements for all entities of dimension `dim'.
  // If `dim' and `tag' are negative, get all the elements in the mesh.
  // `elementTypes' contains the MSH types of the elements (e.g. `2' for 3-node
  // triangles: see `getElementProperties' to obtain the properties for a given
  // element type). `elementTags' is a vector of the same length as
  // `elementTypes'; each entry is a vector containing the tags (unique,
  // strictly positive identifiers) of the elements of the corresponding type.
  // `nodeTags' is also a vector of the same length as `elementTypes'; each
  // entry is a vector of length equal to the number of elements of the given
  // type times the number N of nodes for this type of element, that contains
  // the node tags of all the elements of the given type, concatenated: [e1n1,
  // e1n2, ..., e1nN, e2n1, ...].
  static void get_elements(std::vector<int> & elementTypes,
                           std::vector<std::vector<std::size_t> > & elementTags,
                           std::vector<std::vector<std::size_t> > & nodeTags,
                           const int dim = -1,
                           const int tag = -1);
  // get the number of vertices in a gmsh element type
  static size_t get_n_vertices(const int element_type);
  // get vtk id of a gmsh element type
  static int get_vtk_id(const int element_type);

 private:
  /* Prohibited */
  GmshInterface();
  /* read msh v 2.0 file */
  static void read_msh_v2_(std::fstream & mesh_file, mesh::Mesh & mesh);
  // read .msh v 4.0 file
  static void read_msh_v4_(std::fstream & mesh_file, mesh::Mesh & mesh);
  // build grid for a polyhedron-bounded volume
  static void build_triangulation_(const angem::Polyhedron<double> & cell);
  // estimate the element size for the triangulation element
  static double compute_element_size_(const angem::Polyhedron<double> & cell);
};

}  // end namespace gprs_data
