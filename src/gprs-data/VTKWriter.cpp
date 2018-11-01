#include <VTKWriter.hpp>

#include <fstream>  // stringstream


namespace IO
{

void VTKWriter::write_vtk(const std::vector<Point>                    & vertices,
                          const std::vector<std::vector<std::size_t>> & cells,
                          const std::string                           & fname)
{
  std::ofstream out;
  out.open(fname.c_str());

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

  out.close();
}


void VTKWriter::write_vtk(const std::vector<Point>    & vertices,
                          const std::vector<Gelement> & elements,
                          const std::string           & fname)
{
  std::ofstream out;
  out.open(fname.c_str());

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
  out << std::endl;

  out.close();
}

}  // end namespace
