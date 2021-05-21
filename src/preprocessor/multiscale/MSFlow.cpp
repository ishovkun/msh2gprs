#include "MSFlow.hpp"
#include "MetisInterface.hpp"
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter
#include "EdgeWeightedGraph.hpp"
#include "SupportRegionsFVMGraph.hpp"

namespace multiscale {

using namespace algorithms;

MSFlow::MSFlow(mesh::Mesh const & grid, gprs_data::SimData & data)
    : _grid(grid), _data(data)
{
  size_t const n_coarse = 4;
  PureConnectionMap cell_connections;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face) {
    const auto neighbors = face->neighbors();
    if (neighbors.size() == 2) // not a boundary face
      cell_connections.insert(neighbors[0]->index(), neighbors[1]->index());
  }

  auto partition = MetisInterface<hash_algorithms::empty>::
      build_partitioning(cell_connections, n_coarse, grid.n_active_cells());

  std::string fname = "solution.vtk";
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  mesh::IO::VTKWriter::write_geometry(_grid, out);
  mesh::IO::VTKWriter::enter_section_cell_data(partition.size(), out);
  mesh::IO::VTKWriter::add_data(partition, "partition", out);

  EdgeWeightedGraph g(partition.size());
  for (auto & con: _data.flow_connection_data)
    g.add(UndirectedEdge(con.elements[0], con.elements[1], std::fabs(con.coefficients[0])));

  SupportRegionsFVMGraph regions(partition, std::move(g));
  std::vector<int> centers(partition.size(), 0);
  for ( auto c : regions.centers() )
    centers[c] = 1;

  mesh::IO::VTKWriter::add_data(centers, "centers", out);

  out.close();
}

}  // end namespace multiscale
