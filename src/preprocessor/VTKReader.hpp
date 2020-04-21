#pragma once

#include "mesh/Mesh.hpp"
#include <string>
#include <fstream>  // fstream

namespace io {

/**
 * This class reads a legacy vtk format and assembles
 * mesh::Mesh from it.
 */
class VTKReader
{
 public:
  VTKReader(const std::string file_name, mesh::Mesh & grid);

 protected:
  void read_file_(const std::string & file_name);
  void read_header_(std::fstream & in) const;
  void read_vertices_(std::fstream & in);
  void create_grid_();

  mesh::Mesh & _grid;
  std::vector<angem::Point<3,double>> _vertex_coordinates;
  std::vector<std::vector<size_t>> _cell_entries;
  std::vector<int> _vtk_ids;
};

}  // end namespace gprs_data
