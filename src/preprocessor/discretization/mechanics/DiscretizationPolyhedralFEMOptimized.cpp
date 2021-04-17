#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Elements/PolyhedralElementScaled.hpp"
#include "Isomorphism.hpp"         // provides isomorphism
#include <unordered_map>

namespace discretization {
using Point = angem::Point<3,double>;
using std::vector;

DiscretizationPolyhedralFEMOptimized::
DiscretizationPolyhedralFEMOptimized(mesh::Mesh & grid,
                                     const FiniteElementConfig & config,
                                     const std::vector<int> & fracture_face_markers,
                                     const std::vector<int> & neumann_face_markers)
    : DiscretizationPolyhedralFEM::DiscretizationPolyhedralFEM(grid, config,
                                                               fracture_face_markers,
                                                               neumann_face_markers)

{}

void DiscretizationPolyhedralFEMOptimized::build_(mesh::Cell & cell)
{
  std::vector<size_t> order;  // vertex ordering in case we need to reorder
  if ( auto master = known_element_(cell, order) )
  {
    std::cout << "found similar element " <<  cell.index() << std::endl;
    cell.reorder_vertices(order);
    reorder_faces_(cell, master->host_cell());
    _element = std::make_shared<PolyhedralElementScaled>( cell, _grid, *master, _config );
  }
  else {
    std::cout << "new element topology " <<  cell.index() << std::endl;
    DiscretizationPolyhedralFEM::build_(cell);
    // remember topology for future re-use
    _masters[cell.n_vertices()].push_back(
        std::dynamic_pointer_cast<PolyhedralElementBase>(_element ));
 }
}

std::shared_ptr<PolyhedralElementBase>
DiscretizationPolyhedralFEMOptimized::known_element_(mesh::Cell const & cell,
                                                     std::vector<size_t> & order) const
{
  size_t const hsh = cell.n_vertices();
  if (_masters.count(hsh))
  {
    auto it = _masters.find(hsh);
    for (auto master : it->second)
    {
      auto result = Isomorphism::check(*master->host_topology(),
                                       *cell.polyhedron());
      if (result.first)
      {
        order = result.second;
        return master;
      }
    }
  }
  std::cout << "nada" << std::endl;
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
    order[f] = src_map[ dst_faces[f] ];

  dst.reorder_faces(order);
}

std::vector<uint64_t> DiscretizationPolyhedralFEMOptimized::
compress_faces_(mesh::Cell const & cell) const
{
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

}  // end namespace discretization
