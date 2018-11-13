#include <cell_iterator.hpp>
#include <mesh_methods.hpp>
#include <PolyhedronFactory.hpp>

namespace mesh
{

cell_iterator::
cell_iterator(const std::size_t                       icell,
              angem::PointSet<3,double>             & vertices,
              std::vector<std::vector<std::size_t>> & cells,
              std::unordered_map<uint256_t, std::vector<std::size_t>> & map_faces,
              std::vector<int>                      & shape_ids,
              std::vector<int>                      & cell_markers)
    :
    icell(icell),
    vertices(vertices),
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
  return std::move(get_element_center(vertices, cells[icell]));
}


cell_iterator & cell_iterator::operator++()
{
  icell++;
  return (*this);
}


angem::Polyhedron<double> cell_iterator::polyhedron() const
{
  return angem::PolyhedronFactory::create<double>(vertices.points, cells[icell],
                                                  shape_ids[icell]);
}


}
