#include <VTKWriter.hpp>

#include <fstream>  // stringstream


namespace IO
{
using std::endl;

// edfm
void VTKWriter::write_surface_geometry(const std::vector<Point>    & vertices,
                                       const std::vector<std::vector<std::size_t>> & cells,
                                       std::ofstream               & out)
{
  out << "# vtk DataFile Version 2.0 \n";
  out << "3D Fractures \n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  // points
  const std::size_t n_points = vertices.size();
  out << "POINTS" << "\t"
      << n_points << " float"
      << std::endl;

  for (const auto & p : vertices)
    out << p << std::endl;

  // cells
  const std::size_t n_cells = cells.size();
  std::size_t vind_size_total = 0;
  for (const auto & cell : cells)
    vind_size_total += cell.size();

  out << "CELLS" << "\t"
      << n_cells << "\t"
      << vind_size_total + n_cells
      << std::endl;

  for (const auto & cell : cells)
  {
    out << cell.size() << "\t";
    for (std::size_t i : cell)
      out << i << "\t";
    out << std::endl;
  }

  out << std::endl;
  out << "CELL_TYPES" << "\t" << n_cells << std::endl;
  for (const auto & cell : cells)
  {
    if (cell.size() == 4)  // quad
      out << 9 << std::endl;
    else if (cell.size() == 3)  // triangle
      out << 5 << std::endl;
    else
      out << 7 << std::endl;
    // else if (cell.size() == 5 || cell.size() == 6)  // polygon
    //   out << 7 << std::endl;
    // else
    // {
    //   std::cout << "unknown cell type : " << cell.size() << " vertices" << std::endl;
    //   exit(-1);
    // }
  }

}


void VTKWriter::
write_surface_geometry(const std::vector<Point>                    & vertices,
                       const std::vector<std::vector<std::size_t>> & cells,
                       const std::string                           & fname)
{
  std::ofstream out;
  out.open(fname.c_str());
  write_surface_geometry(vertices, cells, out);
  out.close();
}

void VTKWriter::write_geometry(const Mesh        & grid,
                               const std::string & fname)
{
  std::ofstream out;
  out.open(fname.c_str());
  write_geometry(grid, out);
  out.close();
}


void VTKWriter::write_geometry(const Mesh               & grid,
                               std::ofstream            & out)
{
  if (grid.n_cells() == grid.n_active_cells())
    write_geometry_classic_(grid, out);
  else
    write_geometry_face_based_(grid,out);
}

void VTKWriter::write_geometry_classic_(const Mesh               & grid,
                                        std::ofstream            & out)
{
  out << "# vtk DataFile Version 2.0 \n";
  out << "3D Fractures \n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  // points
  const std::size_t n_points = grid.n_vertices();
  out << "POINTS" << "\t" << n_points << " float" << std::endl;
  for (const auto & p : grid.vertices()) out << p << std::endl;

  // cells
  const std::size_t n_cells = grid.n_active_cells();
  std::size_t vind_size_total = 0;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    vind_size_total += cell->n_vertices();

  out << "CELLS" << "\t"
      << n_cells << "\t"
      << vind_size_total + n_cells
      << std::endl;

  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    out << cell->n_vertices() << "\t";
    for (std::size_t i : cell->vertices())
      out << i << "\t";
    out << std::endl;
  }

  out << std::endl;
  out << "CELL_TYPES" << "\t" << n_cells << std::endl;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    out << cell->vtk_id() << std::endl;
}

void VTKWriter::write_geometry_face_based_(const Mesh & grid,
                                           std::ofstream            & out)
{
  out << "# vtk DataFile Version 2.0" << endl;
  out << "3D Fractures \n";
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID \n";

  // points
  const std::size_t n_points = grid.n_vertices();
  out << "POINTS" << "\t" << n_points << " float" << std::endl;
  for (const auto & p : grid.vertices()) out << p << std::endl;

  const size_t n_entries_total = count_number_of_cell_entries_(grid);

  out << "CELLS " << grid.n_active_cells() << " " << n_entries_total << endl;

  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    out << count_number_of_cell_entries_(*cell) << endl;
    const auto & faces = cell->faces();
    out << faces.size() << endl;
    for (const auto & face : faces)
    {
      const auto & vertices = face->vertices();
      out << vertices.size() << " ";
      for (const size_t v : vertices)
        out << v << " ";
      out << endl;
    }
  }

  out << "CELL_TYPES" << "\t" << grid.n_active_cells() << std::endl;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    out << 42 << endl;
    // out << cell->vtk_id() << std::endl;
}

size_t VTKWriter::count_number_of_cell_entries_(const Mesh & grid)
{
  size_t n = 0;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    n++; // n_entries
    n += count_number_of_cell_entries_(*cell);
  }
  return n;
}

size_t VTKWriter::count_number_of_cell_entries_(const Cell & cell)
{
  size_t n = 0;
  const auto & faces = cell.faces();
  n++; // n_faces
  for (const auto & face : faces)
  {
    n++;  // face.vertices.size()
    n += face->vertices().size();
  }
  return n;
}

// this function is used to write line segments to show wells in paraview
void VTKWriter::write_well_trajectory(const std::vector<Point>                              & vertices,
                                      const std::vector<std::pair<std::size_t,std::size_t>> & indices,
                                      const std::string                                     & fname)
{
  std::ofstream out;
  out.open(fname.c_str());

  out << "# vtk DataFile Version 2.0 \n";
  out << "Wells \n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  // points
  const std::size_t n_points = vertices.size();
  out << "POINTS" << "\t" << n_points << " float" << std::endl;
  for (const auto & p : vertices)
    out << p << std::endl;

  // 3 because n_cells + 2 points per cell
  out << "CELLS" << "\t" << indices.size() << "\t" << 3*indices.size() << std::endl;

  for (const auto & segment : indices)
    out << 2 << "\t" << segment.first << "\t" << segment.second << std::endl;

  out << "CELL_TYPES" << "\t" << indices.size() << std::endl;
  for (std::size_t i=0; i<indices.size(); ++i)
    out << 3 << std::endl;

  out.close();
}


void VTKWriter::enter_section_cell_data(const std::size_t n_cells,
                                        std::ofstream & out)
{
  out << "CELL_DATA" << "\t" << n_cells << std::endl;
}


void VTKWriter::enter_section_point_data(const std::size_t n_vertices,
                                         std::ofstream & out)
{
  out << "POINT_DATA" << "\t" << n_vertices << std::endl;
}

void VTKWriter::write_geometry(const std::vector<angem::Point<3,double>>  &vertices,
                               const std::vector<std::vector<size_t>> &cells,
                               const std::vector<int> & cell_types,
                               std::ofstream & out)
{
  out << "# vtk DataFile Version 2.0 \n";
  out << "3D Element \n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  const std::size_t n_points = vertices.size();
  out << "POINTS" << "\t" << n_points << " float" << std::endl;
  for (const auto & p : vertices) out << p << std::endl;

  // cells
  const std::size_t n_cells = cells.size();
  std::size_t vind_size_total = 0;
  for (const auto & cell : cells)
    vind_size_total += cell.size();

  out << "CELLS" << "\t"
      << n_cells << "\t"
      << vind_size_total + n_cells
      << std::endl;

  for (const auto & cell : cells)
  {
    out << cell.size() << "\t";
    for (const size_t v : cell) out << v << "\t";
    out << std::endl;
  }

  out << endl;
  out << "CELL_TYPES" << "\t" << n_cells << std::endl;
  for (const int type : cell_types)
    out << type << std::endl;
}


}  // end namespace
