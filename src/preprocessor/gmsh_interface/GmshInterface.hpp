#pragma once

#include "element_types.hpp"
#include "mesh/Mesh.hpp"
#include "angem/LineSegment.hpp"
#include <fstream>  // fstream

#include <boost/functional/hash.hpp> //had to add it eventhgough it is bond to gmsh interfaces
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
  static void initialize_gmsh(bool verbose = false);
  // run after wrapping up the work with the grid
  static void finalize_gmsh();
  // build gmsh grid bounded by the faces of the cell
  static void build_triangulation(const mesh::Cell & cell, const double n_vertices_on_edge);
  // build triangulation and save it into grid
  static void build_triangulation(const mesh::Cell & cell, mesh::Mesh & grid,
                                  const double n_vertices_on_edge = 2);
  // build triangulation with embedded points.
  // the outer boundary is specified with angem::Polyhedron.
  // The embedded points are given as a vector of angem::Point
  // The output is written into grid
  static void build_triangulation_embedded_points(angem::Polyhedron<double> const & boundary,
                                                  std::vector<angem::Point<3,double>> const & embedded,
                                                  mesh::Mesh & grid);
  // build triangulation with embedded line segments.
  // the outer boundary is specified with angem::Polyhedron.
  // The embedded lines are given as a vector of angem::LineSegment
  // The output is written into grid
  static void build_triangulation_embedded_lines(angem::Polyhedron<double> const & boundary,
                                                 std::vector<angem::LineSegment<double>> const & embedded,
                                                 mesh::Mesh & grid);
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
  // save current gmsh model
  static void save_gmsh_grid(const std::string fname);
  // get gmsh id from vtk id
  // searches gmsh element types in time linear to the number of gmsh types
  // throws std::invalid_argument
  static int get_gmsh_element_id(const angem::VTK_ID vtk_id);

 private:
  /* Prohibited */
  GmshInterface();
  /* read msh v 2.0 file */
  static void read_msh_v2_(std::fstream & mesh_file, mesh::Mesh & mesh);
  // read .msh v 4.0 file
  static void read_msh_v4_(std::fstream & mesh_file, mesh::Mesh & mesh);
  // extract data from gmsh and save into mesh::Mesh
  // returns vertex numbering
  static std::vector<size_t> extract_grid_(mesh::Mesh & grid);
  // build grid for a polyhedron-bounded volume
  static void build_triangulation_(const angem::Polyhedron<double> & cell,
                                   const double n_vertice_on_edge);
  // estimate the element size for the triangulation element
  static double compute_element_size_(const angem::Polyhedron<double> & cell);
  // conveniance function to insert element into Mesh::mesh
  static void insert_elements_(const int dim, const int tag,
                               const std::vector<size_t> & vertex_numbering,
                               mesh::Mesh & grid);

  static void insert_boundary_data_(angem::Polyhedron<double> const & bnd,
                                    double const n_vertices_on_edge);
  static void insert_vertices_(std::vector<angem::Point<3,double>> const & vertices,
                               std::vector<mesh::Edge> const & edges,
                               double n_vertices_on_edge,
                               int first_tag = 1);
  // convenience function to insert edges into model
  static void insert_boundary_edges_(std::vector<mesh::Edge> const & edges);
  // convenience function to insert faces into model
  static void insert_boundary_faces_(angem::Polyhedron<double> const & bnd);
  // convenience function to insert surfaces into model
  static void insert_surface_loop(size_t n_surfaces);
  static void insert_bounding_volume();
};

}  // end namespace gprs_data
