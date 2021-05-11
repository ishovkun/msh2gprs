#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Elements/PolyhedralElementScaled.hpp"
#include "Isomorphism.hpp"         // provides isomorphism
#include "mesh/io/VTKWriter.hpp"  // debug
#include "logger/Logger.hpp"
#include <unordered_map>
#include <limits>

namespace discretization {
using Point = angem::Point<3,double>;
using std::vector;
static constexpr size_t UNASSIGNED = std::numeric_limits<size_t>::max();

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
  logging::log() << "Enumerating cell types...";
  _element_types.assign( _grid.n_cells_total(), UNASSIGNED );

  std::vector<size_t> order;  // vertex ordering in case we need to reorder
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell) {
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
      _shapes.push_back(std::move(poly));
      _vertex_to_cell_types[cell->n_vertices()].push_back( type );
    }
    _element_types[cell->index()] = type;
  }
  logging::log() << "OK" << "   Various cell types = " << _element_types.size() << std::endl;
}

void DiscretizationPolyhedralFEMOptimized::build_(mesh::Cell & cell)
{
  size_t const type = _element_types[cell.index()];
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

// void DiscretizationPolyhedralFEMOptimized::
// reorder_faces_(mesh::Cell & dst, mesh::Cell const & src) const
// {
//   auto const dst_faces = compress_faces_(dst);
//   auto const src_faces = compress_faces_(src);
//   assert( src_faces.size() == dst_faces.size() && "elements do not match" );

//   size_t const nf = dst_faces.size();
//   std::unordered_map<size_t,size_t> src_map;
//   for (size_t f = 0; f < nf; ++f)
//     src_map[src_faces[f]] = f;

//   std::vector<size_t> order(nf, 0);
//   for (size_t f = 0; f < nf; ++f)
//   {
//     assert( src_map.count(dst_faces[f]) && "isomorphism yielded wrong result; cannot establish face mapping");
//     order[f] = src_map[dst_faces[f]];
//   }

//   auto & faces = dst.face_indices();
//   angem::reorder_to(faces, order);
// }

// std::vector<uint64_t> DiscretizationPolyhedralFEMOptimized::
// compress_faces_(mesh::Cell const & cell) const
// {
//   // Build face-vertex adjacency matrix
//   // instead of using a proper adjacency matrix
//   // we save up on storage by jamming matrix rows
//   // into 64-bit integers
//   assert( cell.n_vertices() < 64 );

//   // vertices: global index to local cell index
//   std::unordered_map<size_t, size_t> global_to_local;
//   auto const verts = cell.vertices();
//   for (size_t v = 0; v < verts.size(); ++v)
//     global_to_local[verts[v]] = v;

//   auto const faces = cell.faces();
//   std::vector<uint64_t> c(faces.size(), 0u);
//   for (size_t f = 0; f < faces.size(); ++f)
//     for (size_t vertex : faces[f]->vertices())
//       c[f] |= (1 << global_to_local[vertex]);

//   return c;
// }

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
