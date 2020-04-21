#include "VTKReader.hpp"

namespace io {

VTKReader::VTKReader(const std::string file_name, mesh::Mesh & grid)
    : _grid(grid)
{
  read_file_(file_name);
}

void VTKReader::read_file_(const std::string & file_name)
{
  std::fstream in;
  in.open(file_name.c_str(), std::fstream::in);
  if (!in) throw std::out_of_range(file_name + " does not exist");

  read_header_(in);
  // read_vertices_(in);
  exit(0);
}

void VTKReader::read_header_(std::fstream & in) const
{
  std::string line;
  std::getline(in, line);  // vtk DataFile Version 3.0
  std::getline(in, line);  // 3D Grid
  std::getline(in, line);  // ASCII
  if (line != "ASCII")
    throw std::invalid_argument("File format not supported");


}

}  // end namespace io
