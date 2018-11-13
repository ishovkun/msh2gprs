#include <mesh_methods.hpp>
#include <angem/utils.hpp>

namespace mesh
{

Point get_element_center(const angem::PointSet<3,double> & vertices,
                         const std::vector<std::size_t>  & ivertices)
{
  return std::move(angem::compute_center_mass
                   (get_vertex_coordinates(vertices, ivertices)));
}


std::vector<Point> get_vertex_coordinates(const angem::PointSet<3,double> & vertices,
                                          const std::vector<std::size_t>  & ivertices)
{
  std::vector<Point> element_vertices(ivertices.size());
  for (std::size_t i=0; i<ivertices.size(); ++i)
    element_vertices[i] = vertices[ivertices[i]];
  return std::move(element_vertices);
}



}
