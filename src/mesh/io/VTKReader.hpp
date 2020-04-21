#pragma once

#include "mesh/Mesh.hpp"
#include <string>
#include <fstream>  // fstream

namespace mesh {

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
  void read_cells_(std::fstream & in);
  void read_cell_types_(std::fstream & in);
  void create_grid_();
  void create_general_polyhedron_cell_(size_t & entry_idx);
  void create_regular_polyhedron_cell_(const int vtk_id, size_t & entry_idx);

  mesh::Mesh & _grid;
  std::vector<size_t> _cell_entries;
  std::vector<int> _vtk_ids;
};

}  // end namespace io

}  // end namespace msh
