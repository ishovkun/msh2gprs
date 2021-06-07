#include "MSFlow.hpp"
#include "MetisInterface.hpp"
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter
#include "algorithms/EdgeWeightedGraph.hpp"
#include "SupportRegionsFVMGraph.hpp"
#include "SupportRegionsLaplaceMod.hpp"
#include "ShapeFunctionSolver.hpp"

namespace multiscale {

using namespace algorithms;

MSFlow::MSFlow(mesh::Mesh const & grid, gprs_data::SimData & data, MultiscaleConfig const &config)
    : _grid(grid)
    , _data(data)
    , _type(config.type)
    ,_ncoarse(config.n_blocks)
{
  size_t const n = data.flow.cv.size();
  auto const weight_func = build_weight_function();
  EdgeWeightedGraph g(n);
  // find min non_zero perm
  for (auto & con: _data.flow.con)
    g.add(UndirectedEdge(con.elements[0], con.elements[1], weight_func(con.coefficients[0])));

  auto partition = MetisInterface::partition(g, _ncoarse);
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
  {
    std::vector<size_t> support(n, 0);
    for (size_t cell : regions.get(coarse))
      support[cell] = 1;
    for (size_t cell : regions.get_boundary(coarse))
      support[cell] = 2;

    for (size_t coarse1 = 0; coarse1 < regions.size(); ++coarse1)
      support[regions.centers()[coarse1]] = 4;

    support[regions.centers()[coarse]] = 3;

    mesh::IO::VTKWriter::add_data(support, "support-" + std::to_string(coarse), out);
  }

  // auto c = regions.centers();
  // for (size_t coarse = 0; coarse < regions.size(); ++coarse) {
  //   std::vector<size_t> bnd = c;
  //   bnd.erase(bnd.begin() + coarse);
  //   ShapeFunctionSolver solver(c[coarse], bnd, _data);
  //   mesh::IO::VTKWriter::add_data(solver.solution(), "support-" + std::to_string(coarse), out);
  // }

  out.close();
}

std::function<double(double)> MSFlow::build_weight_function() const
{
  double minperm = std::numeric_limits<double>::max();
  for (auto & con : _data.flow.con)
    if (con.coefficients[0] != 0.f)
      minperm = std::min(minperm, std::fabs(con.coefficients[0]));
  double const minlogperm = std::log(minperm);
  std::function<double(double)> perm_func = [minlogperm](double k) -> double {
    if (k == 0.f) return 0.f;
    return std::log(std::fabs( k )) - minlogperm + 1;
  };
  return perm_func;
}

}  // end namespace multiscale
