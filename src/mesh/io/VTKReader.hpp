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
  std::vector<std::vector<double>> & get_cell_data() {return _cell_data; }
  std::vector<std::vector<double>> & get_point_data() { return _point_data; }
  std::vector<std::string> & get_cell_data_keys() { return _cell_data_names; }
  std::vector<std::string> & get_point_data_keys() { return _point_data_names; }

 protected:
  void read_file_(const std::string & file_name);
  void read_header_(std::fstream & in) const;
  void read_vertices_(std::fstream & in);
  void read_cells_(std::fstream & in);
  void read_cell_types_(std::fstream & in);
  void create_grid_();
  void create_general_polyhedron_cell_(size_t & entry_idx);
  void create_regular_polyhedron_cell_(const int vtk_id, size_t & entry_idx);
  void read_data_arrays_(std::fstream & in);
  void read_array_(std::fstream & in, std::vector<double> & data,
                   std::string & name, const size_t array_size) const;

  mesh::Mesh & _grid;
  std::vector<size_t> _cell_entries;
  std::vector<int> _vtk_ids;
  std::vector<std::vector<double>> _cell_data;
  std::vector<std::vector<double>> _point_data;
  std::vector<std::string> _cell_data_names;
  std::vector<std::string> _point_data_names;
};

}  // end namespace io

}  // end namespace msh
