#include <cell_iterator.hpp>
#include <mesh_methods.hpp>
#include <PolyhedronFactory.hpp>

namespace mesh
{

cell_iterator::
cell_iterator(const std::size_t                       icell,
              angem::PointSet<3,double>             & vertices,
              std::vector<std::vector<std::size_t>> & cells,
              std::unordered_map<hash_type, Face>   & map_faces,
              std::vector<int>                      & shape_ids,
              std::vector<int>                      & cell_markers)

    :
    icell(icell),
    mesh_vertices(vertices),
    cells(cells),
    map_faces(map_faces),
    shape_ids(shape_ids),
    cell_markers(cell_markers)
{}


bool cell_iterator::operator==(const cell_iterator & other) const
{
  // compare cell index
  if (icell != other.icell)
    return false;
  const std::string msg = "Warning: comparing cell iterators of different mesh objects";

  // compare pointers to data structures
  if ((&map_faces) != (&other.map_faces))
  {
    std::cout << msg << std::endl;
    return false;
  }
  if ((&shape_ids) != (&other.shape_ids))
  {
    std::cout << msg << std::endl;
    return false;
  }
  if ((&cell_markers != (&other.cell_markers)))
  {
    std::cout << msg << std::endl;
    return false;
  }
  return true;
}


bool cell_iterator::operator!=(const cell_iterator & other) const
{
  return !(*this == other);
}


Point cell_iterator::center() const
{
  return get_element_center(mesh_vertices, cells[icell]);
}


cell_iterator & cell_iterator::operator++()
{
  icell++;
  return (*this);
}


std::unique_ptr<Polyhedron<double>> cell_iterator::polyhedron() const
{
  return angem::PolyhedronFactory::create<double>(mesh_vertices.points,
                                                  cells[icell],
                                                  shape_ids[icell]);
}


std::vector<face_iterator> cell_iterator::faces() const
{
  const std::vector<std::vector<std::size_t>> face_vertex_indices =
      angem::PolyhedronFactory::get_global_faces<double>(cells[icell],
                                                         shape_ids[icell]);

  std::vector<face_iterator> faces;
  for (const auto & ivertices : face_vertex_indices)
  {
    const auto hash = hash_value(ivertices);
    auto it_face = map_faces.find(hash);
    if (it_face == map_faces.end())
      throw std::out_of_range("face does not exist");
    face_iterator face(it_face, mesh_vertices);
    faces.push_back(face);
  }
  return faces;
}


bool cell_iterator::has_vertex( const std::size_t ivertex ) const
{
  for (auto & jvertex : cells[icell])
    if (jvertex == ivertex)
      return true;

  return false;
}


std::vector<std::size_t> cell_iterator::neighbor_indices() const
{
  std::vector<std::size_t> neighbors;
  const auto v_faces = faces();
  for (const auto face : v_faces)
    for (const std::size_t neighbor : face.neighbors())
      if (neighbor != icell)
        neighbors.push_back(neighbor);
  return neighbors;
}


double cell_iterator::volume() const
{
  const auto poly =
      angem::PolyhedronFactory::create<double>(mesh_vertices.points,
                                               cells[icell],
                                               shape_ids[icell]);

  return poly->volume();
}

}
