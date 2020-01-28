#pragma once

#include <mesh/Mesh.hpp>
#include <fstream>  // fstream

namespace gprs_data
{

class GmshInterface
{
 public:
  // read .msh input file and save the result into grid
  static void read_msh(const std::string & filename, mesh::Mesh & grid);

  // build gmsh grid bounded by the faces of the cell
  static void build_grid(const mesh::Cell & cell);

 private:
  /* Prohibited */
  GmshInterface();
  /* read msh v 2.0 file */
  static void read_msh_v2_(std::fstream & mesh_file,
                          mesh::Mesh   & mesh);
  // read .msh v 4.0 file
  static void read_msh_v4_(std::fstream & mesh_file,
                                mesh::Mesh   & mesh);
};

}  // end namespace gprs_data
