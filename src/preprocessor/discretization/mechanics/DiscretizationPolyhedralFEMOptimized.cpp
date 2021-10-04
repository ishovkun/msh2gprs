#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Elements/PolyhedralElementScaled.hpp"
#include "Isomorphism.hpp"         // provides isomorphism
#include "mesh/io/VTKWriter.hpp"  // debug
#include "logger/Logger.hpp"
#include <unordered_map>

#ifdef WITH_EIGEN

namespace discretization {
using Point = angem::Point<3,double>;
using std::vector;

DiscretizationPolyhedralFEMOptimized::
DiscretizationPolyhedralFEMOptimized(mesh::Mesh & grid,
                                     const FiniteElementConfig & config,
                                     const std::vector<int> & fracture_face_markers,
                                     const std::vector<size_t> & neumann_face_indices)
    : DiscretizationPolyhedralFEM::DiscretizationPolyhedralFEM(grid, config, fracture_face_markers,
                                                               neumann_face_indices)

{
  _groups.resize( _grid.n_cells_total(), -1 );
}

void DiscretizationPolyhedralFEMOptimized::build_(mesh::Cell & cell)
{
  // std::cout << "\nbuilding cell " << cell.index()
  //           << " (nv = " << cell.n_vertices() << ")" << std::endl;
  std::vector<size_t> order;  // vertex ordering in case we need to reorder
  if ( auto master = known_element_(cell, order) ) {
    // logging::debug() << "found similar element " <<  cell.index() << std::endl;

    auto & verts = cell.vertices();
    angem::reorder_to(verts, order);

    reorder_faces_(cell, master->host_cell());
    _element = std::make_shared<PolyhedralElementScaled>( cell, _grid, *master, _config );
    _groups[ cell.index() ] = _groups[ master->host_cell().index() ];
  }
  else {
    logging::debug() << "new element topology, cell " <<  cell.index() << std::endl;
    DiscretizationPolyhedralFEM::build_(cell);

    // remember topology for future re-use
    _masters[cell.n_vertices()].push_back( std::dynamic_pointer_cast<PolyhedralElementBase>(_element ) );
    _groups[ cell.index() ] = _ngroups++;
 }
}

std::shared_ptr<PolyhedralElementBase>
DiscretizationPolyhedralFEMOptimized::known_element_(mesh::Cell const & cell,
                                                     std::vector<size_t> & order) const
{
  size_t const hsh = cell.n_vertices();

  if (_masters.count(hsh)) {
    auto it = _masters.find(hsh);
    for (auto master : it->second) {
      Isomorphism isomorphism(*master->host_topology(), *cell.polyhedron());
      if (isomorphism.check())
      {
        order = isomorphism.ordering();
        return master;
      }
    }
  }

  return nullptr;
}

void DiscretizationPolyhedralFEMOptimized::
reorder_faces_(mesh::Cell & dst, mesh::Cell const & src) const
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
    if ( src_map.count(dst_faces[f]) == 0 )
    {
    std::cout << "src faces: ";
    for (auto c : src_faces) std::cout << c << " ";
    std::cout << std::endl;
    std::cout << "dst faces: ";
    for (auto c : dst_faces) std::cout << c << " ";
    std::cout << std::endl;
      std::cout << "src map" << std::endl;
      for (auto [c, ff] : src_map)
        std::cout << ff << " " << c << std::endl;
      std::cout << "dst face " << dst_faces[f] << std::endl;
      mesh::IO::VTKWriter::write_geometry(_grid, dst, "output/dst-" + std::to_string(dst.index()) + ".vtk");
      mesh::IO::VTKWriter::write_geometry(_grid, src, "output/src-" + std::to_string(src.index()) + ".vtk");
    }
    assert( src_map.count(dst_faces[f]) );
    order[f] = src_map[dst_faces[f]];
  }

  auto & faces = dst.face_indices();
  angem::reorder_to(faces, order);
}

std::vector<uint64_t> DiscretizationPolyhedralFEMOptimized::
compress_faces_(mesh::Cell const & cell) const
{
  // Build face-vertex adjacency matrix
  // instead of using a proper adjacency matrix
  // we save up on storage by jamming matrix rows
  // into 64-bit integers
  assert( cell.n_vertices() < 64 );

  // vertices: global index to local cell index
  std::unordered_map<size_t, size_t> global_to_local;
  auto const verts = cell.vertices();
  for (size_t v = 0; v < verts.size(); ++v)
    global_to_local[verts[v]] = v;

  auto const faces = cell.faces();
  std::vector<uint64_t> c(faces.size(), 0u);
  for (size_t f = 0; f < faces.size(); ++f)
    for (size_t vertex : faces[f]->vertices())
      c[f] |= (1 << global_to_local[vertex]);

  return c;
}

std::vector<int> DiscretizationPolyhedralFEMOptimized::get_cell_isomorphic_groups() const
{
  // std::vector<int> groups( _grid.n_cells_total() );
  // int group = 0;
  // std::vector<size_t> order;  // vertex ordering in case we need to reorder
  // for (auto [cell, counter] = std::tuple( _grid.begin_active_cells(), 0 ); cell != _grid.end_active_cells(); ++cell, ++counter) {
  //   if ( auto master = known_element_(*cell, order) ) {
  //     groups[ cell->index() ] = groups[ master->host_cell().index() ];
  //   }
  //   else groups[ cell->index() ] = group++;
  // }

    std::vector<int> ans( _grid.n_active_cells() );
    for (auto [cell, counter] = std::tuple( _grid.begin_active_cells(), 0 ); cell != _grid.end_active_cells(); ++cell, ++counter)
      ans[counter] = _groups[ cell->index() ];
    return ans;
}

}  // end namespace discretization

#endif
