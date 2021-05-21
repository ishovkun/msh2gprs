#include "Idea.hpp"
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter
#include "ShapeFunctionSolver.hpp"
#include <string>
#include <fstream>
#include <algorithm>  // std::swap

namespace multiscale {


Idea::Idea(mesh::Mesh const & grid, gprs_data::SimData & data)
    : _grid(grid), _data(data)
{
  // // size_t source = find_center_cell_();
  // size_t source = 20;
  // // auto bnd = find_boundary_cells_();
  // auto bnd = find_boundary_cells_({ 2, 3 });
  // std::cout << "bnd: ";
  // for (auto c : bnd)
  //   std::cout << c << " ";
  // std::cout << std::endl;
  // ShapeFunctionSolver solver(source, bnd, _grid, _data);
  // auto solution = solver.solution();
  // debug_save_solution_("solution.vtk", solution);

  std::vector<double> sol1, sol2, sol3, sol4;
  {  // top left corner source
    std::cout << "building 1" << std::endl;
    size_t source = 2550;
    auto bnd = find_boundary_cells_({3, 2});
    ShapeFunctionSolver solver(source, bnd, _grid, _data);
    sol1 = solver.solution();
    // debug_save_solution_("solution0.vtk", solver.solution());
  }

  {  // top right corner source
    std::cout << "building 2" << std::endl;
    size_t source = 2600;
    // size_t source = 24;
    auto bnd = find_boundary_cells_({3, 1});
    ShapeFunctionSolver solver(source, bnd, _grid, _data);
    // debug_save_solution_("solution1.vtk", solver.solution());
    sol2 = solver.solution();
  }

  {  // bottom left
    std::cout << "building 3" << std::endl;
    size_t source = 0;
    auto bnd = find_boundary_cells_({4, 2});
    ShapeFunctionSolver solver(source, bnd, _grid, _data);
    // debug_save_solution_("solution2.vtk", solver.solution());
    sol3 = solver.solution();
  }

  {  // bottom left
    std::cout << "building 4" << std::endl;
    size_t source = 50;
    // size_t source = 4;
    auto bnd = find_boundary_cells_({4, 1});
    ShapeFunctionSolver solver(source, bnd, _grid, _data);
    // debug_save_solution_("solution3.vtk", solver.solution());
    sol4 = solver.solution();
  }
  std::string fname = "solution.vtk";
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  mesh::IO::VTKWriter::write_geometry(_grid, out);
  mesh::IO::VTKWriter::enter_section_cell_data(_grid.n_active_cells(), out);
  mesh::IO::VTKWriter::add_data(sol1, "solution1", out);
  mesh::IO::VTKWriter::add_data(sol2, "solution2", out);
  mesh::IO::VTKWriter::add_data(sol3, "solution3", out);
  mesh::IO::VTKWriter::add_data(sol4, "solution4", out);
  std::vector<double> sum(sol1.size(), 0);
  for (size_t i = 0; i < sol1.size(); ++i)
    sum[i] = sol1[i] + sol2[i] + sol3[i] + sol4[i];
  mesh::IO::VTKWriter::add_data(sum, "sum", out);

  std::vector<double> sol1n(sol1.size(), 0);
  std::vector<double> sol2n(sol1.size(), 0);
  std::vector<double> sol3n(sol1.size(), 0);
  std::vector<double> sol4n(sol1.size(), 0);
  for (size_t i = 0; i < sol1.size(); ++i)
  {
    sol1n[i] = sol1[i] / sum[i];
    sol2n[i] = sol2[i] / sum[i];
    sol3n[i] = sol3[i] / sum[i];
    sol4n[i] = sol4[i] / sum[i];
  }
  mesh::IO::VTKWriter::add_data(sol1n, "solution1n", out);
  mesh::IO::VTKWriter::add_data(sol2n, "solution2n", out);
  mesh::IO::VTKWriter::add_data(sol3n, "solution3n", out);
  mesh::IO::VTKWriter::add_data(sol4n, "solution4n", out);
  out.close();

}

std::vector<size_t> Idea::find_boundary_cells_(std::vector<int> const & markers )
{
  std::unordered_set<size_t> bnd;
  for (auto face = _grid.begin_active_faces(); face != _grid.end_active_faces(); ++face)
    if ( std::find(markers.begin(), markers.end(), face->marker()) != markers.end() )
      bnd.insert(face->neighbors()[0]->index());
  return std::vector<size_t>(bnd.begin(), bnd.end());
}

void Idea::debug_save_solution_(std::string const & fname, std::vector<double> const& vec) const
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  mesh::IO::VTKWriter::write_geometry(_grid, out);
  mesh::IO::VTKWriter::enter_section_cell_data(_grid.n_active_cells(), out);
  mesh::IO::VTKWriter::add_data(vec, "solution", out);
  out.close();
}

size_t Idea::find_center_cell_() const
{
  angem::Point<3,double> c;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    c += cell->center();
  c /= _grid.n_active_cells();

  size_t s = 0;
  double mindist = std::numeric_limits<double>::max();
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    double const dist = cell->center().distance(c);
    if (dist < mindist)
    {
      mindist = dist;
      s = cell->index();
    }
  }
  return s;
}

std::vector<size_t> Idea::find_boundary_cells_() const
{
  std::vector<size_t> ans;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    if ( cell->neighbors().size() < 4 )  // cartesian only
      ans.push_back(cell->index());
  }
  return ans;
}

}  // end namespace multiscale
