#include "Subdivision.hpp"

namespace mesh {

Subdivision::Subdivision(const mesh::Cell & cell, mesh::Mesh & triangulation, const size_t order)
    : _parent_cell(cell), _grid(triangulation), _order(order)
{
  if (!_grid.empty())
    throw std::invalid_argument("Grid should be empty");

  // copy cell from its master grid to a new grid
  // this cell will be the parent of all cells in the grid
  create_master_cell_();

  perform_subdivision_(*_grid.begin_active_cells());

  for (size_t i=0; i<order; ++i)  // loop subdivision orders
  {
    // save current cell indices
    std::vector<size_t> active_cells(_grid.n_active_cells());
    size_t j = 0;
    for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
      active_cells[j] = cell->index();
    for (const size_t icell : active_cells)
      perform_subdivision_(_grid.cell(icell));
  }
}

void Subdivision::create_master_cell_()
{
  _grid.vertices() = _parent_cell.vertex_coordinates();
  std::map<size_t,size_t> old_to_new;
  size_t i = 0;
  for (const size_t v : _parent_cell.vertices())
    old_to_new[v] = i++;
  const auto faces = _parent_cell.faces();
  std::vector<FaceTmpData> new_faces(faces.size());
  i = 0;
  for (const auto face : faces)
  {
    const std::vector<std::size_t> & vertices = face->vertices();
    auto & f = new_faces[i];
    f.vertices.reserve(vertices.size());
    for (size_t v=0; v<vertices.size(); ++v)
      f.vertices.push_back( old_to_new[vertices[v]] );
    f.vtk_id = face->vtk_id();
    f.marker = i + 1;  // to comply with old face markering by gmsh
    i++;
  }

  std::vector<size_t> use_faces(new_faces.size());
  std::iota( use_faces.begin(), use_faces.end(), 0 );
  std::vector<size_t> ivertices(_grid.n_vertices());
  std::iota( ivertices.begin(), ivertices.end(), 0 );
  std::cout << std::flush << std::endl;
  _grid.insert_cell_(ivertices, use_faces, new_faces,
                     _parent_cell.vtk_id(), _parent_cell.marker());
}

void Subdivision::perform_subdivision_(mesh::Cell & cell)
{
  const size_t cell_center_index = _grid.n_vertices();
  _grid.vertices().push_back( cell.center() );

  // save face indices since inserting new cells invalidates pointsers
  std::vector<size_t> face_indices;
  for( auto face :cell.faces() )
    face_indices.push_back( face->index() );

  for (const size_t iface : face_indices)
  {
    const mesh::Face * face = &_grid.face(iface);
    const auto c = face->center();
    size_t face_center_index;
    const size_t v = _created_vertices.find(c);
    if ( v == _created_vertices.size() )
    {
      _created_vertices.insert(c);
      _created_vertex_indices.push_back(_created_vertices.size() - 1);
      face_center_index = _grid.insert_vertex(c);
    }
    else face_center_index = _created_vertex_indices[v];

    auto build_trgl_face = [](const std::vector<size_t> vertices, const size_t parent,
                              const int marker, FaceTmpData & f)
                           {
                             f.vertices = vertices;
                             f.vtk_id = angem::VTK_ID::TriangleID;
                             f.parent = parent;
                             f.marker = marker;
                           };

    // refine face and build tetras
    const auto face_vertices = face->vertices();
    for (size_t i=0; i<face_vertices.size(); ++i)
    {
      std::cout << "\n\n " << std::endl;
      size_t i1 = face_vertices[i], i2;
      if (i == face_vertices.size() - 1) i2 = face_vertices[0];
      else                               i2 = face_vertices[i+1];
      std::vector<FaceTmpData> tetra_faces(4);
      // base of the tetra resides on the parent face
      build_trgl_face({i1, i2, face_center_index}, face->index(), face->marker(), tetra_faces[0]);
      // add three more faces to build the tetrahedron
      build_trgl_face({i1, i2, cell_center_index}, constants::invalid_index,
                      constants::default_face_marker, tetra_faces[1]);
      build_trgl_face({i1, face_center_index, cell_center_index}, constants::invalid_index,
                      constants::default_face_marker, tetra_faces[2]);
      build_trgl_face({i2, face_center_index, cell_center_index}, constants::invalid_index,
                      constants::default_face_marker, tetra_faces[2]);
      std::vector<size_t> take_faces(4);
      std::iota(take_faces.begin(), take_faces.end(), 0);
      const size_t child_cell_index =
          _grid.insert_cell_({i1, i2, face_center_index, cell_center_index},
                             take_faces, tetra_faces, angem::TetrahedronID,
                             cell.marker());
      std::cout << "child_cell_index = " << child_cell_index << std::endl;
      cell.m_children.push_back(child_cell_index);
      _grid.m_cells[child_cell_index].m_parent = cell.index();
    }
  }
}

}  // end namespace discretization
