#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Elements/PolyhedralElementScaled.hpp"
#include "Isomorphism.hpp"         // provides isomorphism
#include "mesh/io/VTKWriter.hpp"  // debug
#include "logger/Logger.hpp"
#include "logger/ProgressBar.hpp"  // provides ProgressBar
#include "GlobalOpts.hpp"
#include <unordered_map>
#include <limits>

namespace discretization {
using Point = angem::Point<3,double>;
using std::vector;
static constexpr int UNASSIGNED = -1;

DiscretizationPolyhedralFEMOptimized::
DiscretizationPolyhedralFEMOptimized(mesh::Mesh & grid,
                                     const FiniteElementConfig & config,
                                     const std::vector<int> & fracture_face_markers,
                                     const std::vector<size_t> & neumann_face_indices)
    : DiscretizationPolyhedralFEM::DiscretizationPolyhedralFEM(grid, config, fracture_face_markers,
                                                               neumann_face_indices)

{
  enumerate_elements_();
  _master_elements.assign( _shapes.size(), nullptr);
}

void DiscretizationPolyhedralFEMOptimized::enumerate_elements_()
{
  std::unique_ptr<logging::ProgressBar> progress = nullptr;
  if (gprs_data::GlobalOpts::ref().print_progressbar)
    progress = std::make_unique<logging::ProgressBar>("Enumerating cell types", _grid.n_active_cells());
  else logging::log() << "Enumerating cell types...";

  _cell_types.assign( _grid.n_cells_total(), UNASSIGNED );

  std::vector<size_t> order;  // vertex ordering in case we need to reorder
  for (auto [cell, counter] = std::tuple( _grid.begin_active_cells(), 0 ); cell != _grid.end_active_cells(); ++cell, ++counter) {
    if (progress) progress->set_progress(counter);

    auto poly = cell->polyhedron();
    size_t const type = known_element_(*poly, order);

    if ( type < _shapes.size() ) {
      // reorder cell vertices in order of master's vertices
      angem::reorder_to(cell->vertices(), order);
      // reorder cell faces in order of master's vertices
      // NOTE: here we call cell->polyhedron again, since vertex nubmering changed
      auto face_order = get_face_ordering_(*cell->polyhedron(), *_shapes[type]);
      angem::reorder_to(cell->face_indices(), face_order);
    }
    else {  // new element
      // std::cout << "cell " << cell->index() << " is new" << std::endl;
      _shapes.push_back(std::move(poly));
      _vertex_to_cell_types[cell->n_vertices()].push_back( type );
    }
    _cell_types[cell->index()] = type;
  }

  if (progress) progress->finalize();
  else logging::log() << "OK." << std::endl;

  logging::log() << " Unique cell types = " << _shapes.size() << std::endl;
}

void DiscretizationPolyhedralFEMOptimized::build_(mesh::Cell & cell)
{
  size_t const type = _cell_types[cell.index()];
  if (_master_elements[type])
    _element = std::make_shared<PolyhedralElementScaled>( cell, _grid, *_master_elements[type], _config );
  else {
    DiscretizationPolyhedralFEM::build_(cell);
    _master_elements[type] = std::dynamic_pointer_cast<PolyhedralElementBase>(_element );
  }
}

size_t DiscretizationPolyhedralFEMOptimized::known_element_(angem::Polyhedron<double> const & poly, std::vector<size_t> & order) const
{
  size_t const hash = poly.get_points().size();  // hash

  if ( _vertex_to_cell_types.count(hash) ) {
    auto const it = _vertex_to_cell_types.find(hash);
    for (size_t const type : it->second) {
      Isomorphism isomorphism(*_shapes[type], poly);
      if (isomorphism.check())
      {
        order = isomorphism.ordering();
        return type;
      }
    }
  }

  return _shapes.size();
}

std::vector<size_t> DiscretizationPolyhedralFEMOptimized::
get_face_ordering_(angem::Polyhedron<double> const & dst, angem::Polyhedron<double> const & src) const
{
  auto const dst_faces = compress_faces_(dst);
  auto const src_faces = compress_faces_(src);
  assert( src_faces.size() == dst_faces.size() && "elements do not match" );

  size_t const nf = dst_faces.size();
  std::unordered_map<size_t,size_t> src_map;
  for (size_t f = 0; f < nf; ++f)
    src_map[src_faces[f]] = f;

  std::vector<size_t> order(nf, 0);
  for (size_t f = 0; f < nf; ++f)
  {
    assert( src_map.count(dst_faces[f]) && "isomorphism yielded wrong result; cannot establish face mapping");
    order[f] = src_map[dst_faces[f]];
  }
  return order;
}

std::vector<uint64_t> DiscretizationPolyhedralFEMOptimized::
compress_faces_(angem::Polyhedron<double> const & cell) const
{
  // Build face-vertex adjacency matrix
  // instead of using a proper adjacency matrix

  /* Instead of building a real matrix, we save up on memory and comparison operations
   * by jamming matrix rows into 64-bit integers
   * Example:
   * | 0 1 1 |      |3| (binary 011)
   * | 1 0 1 | ---> |5| (binary 101)
   * | 0 0 1 |      |1| (binary 001)
  */
  auto const & vertices = cell.get_points();
  assert( vertices.size() < 64 );

  auto const & faces = cell.get_faces();
  std::vector<uint64_t> compressed(faces.size(), 0u);
  for (size_t f = 0; f < faces.size(); ++f)
    for (size_t const vertex : faces[f])
      compressed[f] |= (1 << vertex);

  return compressed;
}

}  // end namespace discretization
