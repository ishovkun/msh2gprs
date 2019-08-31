#include <VTKWriter.hpp>

#include <fstream>  // stringstream


namespace IO
{

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

void VTKWriter::write_geometry(const std::vector<Point>    & vertices,
                               const std::vector<Gelement> & elements,
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
  const std::size_t n_cells = elements.size();
  std::size_t vind_size_total = 0;
  for (const auto & cell : elements)
    vind_size_total += cell.vVertices.size();

  out << "CELLS" << "\t"
      << n_cells << "\t"
      << vind_size_total + n_cells
      << std::endl;

  for (const auto & cell : elements)
  {
    out << cell.vVertices.size() << "\t";

    if(cell.vtkIndex == 25) // super wierd element 25
    {
      for (int j = 0; j < 8; j++)
        out << cell.vVertices[j] << "\t";

      out << cell.vVertices[8]  << "\t";
      out << cell.vVertices[11] << "\t";
      out << cell.vVertices[13] << "\t";
      out << cell.vVertices[9]  << "\t";

      out << cell.vVertices[16] << "\t";
      out << cell.vVertices[18] << "\t";
      out << cell.vVertices[19] << "\t";
      out << cell.vVertices[17] << "\t";

      out << cell.vVertices[10] << "\t";
      out << cell.vVertices[12] << "\t";
      out << cell.vVertices[14] << "\t";
      out << cell.vVertices[15] << "\t";
      out << std::endl;
    }
    else if (cell.vtkIndex == 26)
    {
      for (int j = 0; j < 6; j++)
        out << cell.vVertices[j] << "\t";

      out << cell.vVertices[6] << "\t";
      out << cell.vVertices[9] << "\t";
      out << cell.vVertices[7] << "\t";

      out << cell.vVertices[12] << "\t";
      out << cell.vVertices[14] << "\t";
      out << cell.vVertices[13] << "\t";

      out << cell.vVertices[8]  << "\t";
      out << cell.vVertices[10] << "\t";
      out << cell.vVertices[11] << "\t";
      out << std::endl;
    }
    else
    {
      for (int j = 0; j < cell.vVertices.size(); j++)
        out << cell.vVertices[j] << "\t";
      out << std::endl;
    }

  }

  // cell type
  out << "CELL_TYPES" << "\t" << n_cells << std::endl;
  for (const auto & cell : elements)
    out <<  cell.vtkIndex << std::endl;

}


void VTKWriter::write_geometry(const Mesh               & grid,
                               const std::string               & fname)
{
  std::ofstream out;
  out.open(fname.c_str());
  write_geometry(grid, out);
  out.close();
}


void VTKWriter::write_geometry(const Mesh               & grid,
                               std::ofstream            & out)
{
  out << "# vtk DataFile Version 2.0 \n";
  out << "3D Fractures \n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  // points
  const std::size_t n_points = grid.n_vertices();
  out << "POINTS" << "\t"
      << n_points << " float"
      << std::endl;

  for (const auto & p : grid.get_vertices())
    out << p << std::endl;

  // cells
  const std::size_t n_cells = grid.n_cells();
  std::size_t vind_size_total = 0;
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    vind_size_total += cell.n_vertices();

  out << "CELLS" << "\t"
      << n_cells << "\t"
      << vind_size_total + n_cells
      << std::endl;

  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    out << cell.n_vertices() << "\t";
    for (std::size_t i : cell.vertex_indices())
      out << i << "\t";
    out << std::endl;
  }

  out << std::endl;
  out << "CELL_TYPES" << "\t" << n_cells << std::endl;
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    out << cell.vtk_id() << std::endl;
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

}  // end namespace
