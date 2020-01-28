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
