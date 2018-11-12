#pragma once

#include <mesh/Mesh.hpp>

namespace Parsers
{


class GmshReader
{
 public:
  static void read_input(const std::string & filename,
                         mesh::Mesh        & mesh);

  typedef std::unordered_map<int,int> MapIntInt;
  static int get_vtk_index(const int gmsh_index) {return map_gmsh_vtk[gmsh_index];}


 private:
  GmshReader();
  static void read_gmsh2_input(std::fstream & mesh_file,
                               mesh::Mesh   & mesh);
  // variables
  static MapIntInt map_gmsh_vtk;
};


}
