#include "MSFlow.hpp"
#include "MetisInterface.hpp"
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter
#include "algorithms/EdgeWeightedGraph.hpp"
#include "SupportRegionsFVMGraph.hpp"
#include "SupportRegionsLaplaceMod.hpp"
#include "ShapeFunctionSolver.hpp"

namespace multiscale {

using namespace algorithms;

MSFlow::MSFlow(mesh::Mesh const & grid, gprs_data::SimData & data)
    : _grid(grid), _data(data)
{
  size_t const n_coarse = 16;

  EdgeWeightedGraph g(data.cv_data.size());
  for (auto & con: _data.flow_connection_data)
    g.add(UndirectedEdge(con.elements[0], con.elements[1], std::fabs(con.coefficients[0])));

  auto partition = MetisInterface::partition(g, n_coarse);
  std::string fname = "solution.vtk";
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  mesh::IO::VTKWriter::write_geometry(_grid, out);
  mesh::IO::VTKWriter::enter_section_cell_data(partition.size(), out);
  mesh::IO::VTKWriter::add_data(partition, "partition", out);

  SupportRegionsFVMGraph regions(partition, std::move(g));
  // SupportRegionsLaplaceMod regions(partition, data);
  std::vector<int> centers(partition.size(), 0);
  for ( auto c : regions.centers() )
    centers[c] = 1;

  mesh::IO::VTKWriter::add_data(centers, "centers", out);
  for (size_t coarse = 0; coarse < regions.size(); ++coarse)
    mesh::IO::VTKWriter::add_data(regions.get(coarse), "support-" + std::to_string(coarse), out);

  // auto c = regions.centers();
  // for (size_t coarse = 0; coarse < regions.size(); ++coarse) {
  //   std::vector<size_t> bnd = c;
  //   bnd.erase(bnd.begin() + coarse);
  //   ShapeFunctionSolver solver(c[coarse], bnd, _data);
  //   mesh::IO::VTKWriter::add_data(solver.solution(), "support-" + std::to_string(coarse), out);
  // }

  out.close();
}

}  // end namespace multiscale
