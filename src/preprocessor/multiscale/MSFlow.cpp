#include "MSFlow.hpp"
#include "MetisInterface.hpp"
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter
#include "algorithms/EdgeWeightedGraph.hpp"
#include "SupportRegionsFVMGraph.hpp"
#include "SupportRegionsLaplaceMod.hpp"
#include "ShapeFunctionSolver.hpp"
#include "GeometricPartition.hpp"

namespace multiscale {

using namespace algorithms;

MSFlow::MSFlow(mesh::Mesh const & grid, gprs_data::SimData & data, MultiscaleConfig const &config)
    : _grid(grid)
    , _data(data)
    // , _type(config.support_type)
    , _ncoarse( std::accumulate(config.n_blocks.begin(), config.n_blocks.end(), 1, std::multiplies<size_t>()) )
{
  size_t const n = data.flow.cv.size();
  auto const weight_func = build_weight_function();
  EdgeWeightedGraph g(n);
  // find min non_zero perm
  for (auto & con: _data.flow.con)
    g.add(UndirectedEdge(con.elements[0], con.elements[1], weight_func(con.coefficients[0])));

  // partition
  if ( config.part_type == MSPartitioning::metis )
    _part = MetisInterface::partition(g, _ncoarse);
  else if (config.part_type == MSPartitioning::geometric)
    _part = GeometricPartition(config.n_blocks, data).get();

  std::string fname = "solution.vtk";
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  mesh::IO::VTKWriter::write_geometry(_grid, out);
  mesh::IO::VTKWriter::enter_section_cell_data(_part.size(), out);
  mesh::IO::VTKWriter::add_data(_part, "partition", out);

  SupportRegionsFVMGraph regions(_part, std::move(g));
  // SupportRegionsLaplaceMod regions(partition, data);

  _centers = regions.centers();
  // std::vector<int> centers(_part.size(), 0);
  // for ( auto c : regions.centers() )
  //   centers[c] = 1;
  // mesh::IO::VTKWriter::add_data(centers, "centers", out);

  _support.resize( _ncoarse );
  _bnds.resize( _ncoarse );
  for (size_t coarse = 0; coarse < regions.size(); ++coarse)
  {
    _support[coarse] = regions.get(coarse);
    _bnds[coarse] = regions.get_boundary(coarse);

    // visualization output
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

void MSFlow::fill_output_model(MultiScaleOutputData & model) const
{
  model.cell_data = true;
  model.n_coarse = _ncoarse;
  // partition
  model.partitioning.resize(_part.size());
  copy(_part.cbegin(), _part.cend(), model.partitioning.begin());
  // centroids
  model.centroids = _centers;
  // internal
  model.support_internal = _support;
  model.support_boundary = _bnds;
}


}  // end namespace multiscale
