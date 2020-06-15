#pragma once

#include "angem/Point.hpp"
#include "Mesh.hpp"
#include <fstream>

namespace IO
{

using Point = angem::Point<3,double>;
using Mesh = mesh::Mesh;
using mesh::Cell;

class VTKWriter
{
 public:
  // for fractures
  static void write_surface_geometry(const std::vector<Point>                    & vertices,
                                     const std::vector<std::vector<std::size_t>> & cells,
                                     const std::string                           & fname);
  // for fractures
  static void write_surface_geometry(const std::vector<Point>                    & vertices,
                                     const std::vector<std::vector<std::size_t>> & cells,
                                     std::ofstream                               & out);
  // for gmsh
  static void write_geometry(const std::vector<angem::Point<3,double>>  &vertices,
                             const std::vector<std::vector<size_t>> &cells,
                             const std::vector<int> & cell_types,
                             std::ofstream & out);

  static void write_geometry(const Mesh & grid, const std::string & fname);
  static void write_geometry(const Mesh & grid, std::ofstream & out);

  // save geometry of a single cell
  static void write_geometry(const Mesh             & grid,
                             const Cell             & cell,
                             std::string            file_name);

  static void write_well_trajectory(const std::vector<Point>                              & vertices,
                                    const std::vector<std::pair<std::size_t,std::size_t>> & indices,
                                    const std::string                                     & fname);
  // add cell data to vtk file
  template <typename T>
  static void add_data(const std::vector<T> & property,
                       const std::string           keyword,
                       std::ofstream             & out);

  static void enter_section_cell_data(const std::size_t n_cells,
                                      std::ofstream & out);
  static void enter_section_point_data(const std::size_t n_vertices,
                                       std::ofstream & out);
  static size_t count_number_of_cell_entries_(const Mesh & grid);
  static size_t count_number_of_cell_entries_(const Cell & cell);

  template <typename T>
  static void write(const angem::Polyhedron<T> & polyhedron, const std::string & fname);

 protected:
  static void write_geometry_classic_(const Mesh & grid, std::ofstream & out);
  static void write_geometry_face_based_(const Mesh & grid, std::ofstream & out);

 private:
  VTKWriter();
};

// add cell data to vtk file
template <typename T>
void VTKWriter::add_data(const std::vector<T> &     property,
                         const std::string          keyword,
                         std::ofstream            & out)
{
  out << "SCALARS\t" << keyword << "\t";
  out << "float" << std::endl;
  out << "LOOKUP_TABLE HSV" << std::endl;
  for (const double item : property)
    out << static_cast<double>(item)<< std::endl;
  out << std::endl;
}

template <typename T>
void VTKWriter::write(const angem::Polyhedron<T> & polyhedron, const std::string & fname)
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  out << "# vtk DataFile Version 2.0 \n";
  out << "3D Polyhedron \n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  const std::size_t n_points = polyhedron.get_points().size();
  out << "POINTS" << "\t" << n_points << " float" << "\n";
  for (const auto & p : polyhedron.get_points()) out << p << "\n";

  size_t n_entries = 1;  // << n_faces
  for (const auto & face : polyhedron.get_faces())
    n_entries += face.size() + 1;  // << n_face_vertices + each face vertex
  out << "CELLS " << 1 << " " << n_entries + 1 << "\n";

  out << n_entries << "\n";
  out << polyhedron.get_faces().size() << "\n";
  for (const auto & face : polyhedron.get_faces())
  {
    out << face.size() << " ";
    for (const size_t v : face)
      out << v << " ";
    out << "\n";
  }
  out << "CELL_TYPES" << "\t" << 1 << "\n";
  out << angem::GeneralPolyhedronID << "\n";

  out.close();
}

}  // end namespace
