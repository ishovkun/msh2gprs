#include "CellSplitter.hpp"
#include "angem/Collisions.hpp"    // angem::split

namespace mesh {

CellSplitter::CellSplitter(Mesh & grid)
    : _grid(grid)
{}

void CellSplitter::split_cell(Cell cell, const angem::Plane<double> & plane,
                              const int splitting_face_marker)
{
  if (!_grid.cell(cell.index()).is_active())
  {
    const auto children = _grid.cell(cell.index()).immediate_children();
    // make sure this is a cell with hanging nodes
    assert ( children.size() == 1 );
    return split_cell(*children[0], plane, splitting_face_marker);
  }
  // std::cout << splitting_face_marker<< "-split " << cell.index() << " (parent "
  //           << cell.m_parent << " ult " << cell.ultimate_parent().index() << ")"
  //           << std::endl << std::flush;
  assert (cell.is_active());
  const size_t parent_cell_index = cell.index();

  // Bookkeeping:
  //  fill polygroup's internal set with the existing vertex coordinates
  // in order to have a map of those to the global vertex indices,
  // which will come in handy when inserting new splitted cells into grid.
  // We can do it because splitting will insert the same vertices plus
  // those that appeared due to plase-face intersection.
  angem::PolyGroup<double> split;
  std::vector<size_t> global_vertex_indices;
  for (const size_t vertex : cell.vertices())
  {
    split.vertices.insert(_grid.vertex(vertex));
    global_vertex_indices.push_back(vertex);
  }
  // the actual geometry happens here
  const std::unique_ptr<angem::Polyhedron<double>> polyhedron = cell.polyhedron();
  const double h = cell.center().distance(cell.vertex_coordinates().front());

  // which polygons in split belong to which faces
  std::vector<size_t> polygroup_polygon_parents;
  const bool success = angem::split(*polyhedron, plane, split, polygroup_polygon_parents,
               constants::marker_below_splitting_plane,
               constants::marker_above_splitting_plane,
               constants::marker_splitting_plane,
               /* tol = */ 1e-4*h);
  if (!success) return;

  const std::vector<Face*> & cell_faces = cell.faces();

  // insert new vertices (those that occured due to splitting)
  std::vector<size_t> new_vertices;
  track_new_vertices_(global_vertex_indices, split, new_vertices);

  // map local indices to global
  std::vector<std::vector<size_t>> face_vertex_global_numbering;
  for (size_t i = 0; i < split.polygons.size(); i++)
    face_vertex_global_numbering.push_back(build_global_face_indices_
                                           (split.polygons[i], global_vertex_indices));

  // find the index of the splitting face
  size_t split_face_local_index;  // need this to figure out parent/child faces
  for (size_t i = 0; i < split.polygons.size(); i++)
    if ( split.markers[i] == constants::marker_splitting_plane )
      split_face_local_index = i;

  // make two groups of faces (polygons that constitute polyhedra) that will form the new cells
  std::vector<FaceTmpData> tmp_faces(split.polygons.size());
  std::vector<size_t> cell_above_faces, cell_below_faces;
  create_face_groups_(split, face_vertex_global_numbering, polygroup_polygon_parents,
                      cell_faces, splitting_face_marker, split_face_local_index,
                      cell_above_faces, cell_below_faces, tmp_faces);

  // insert new cells
  const size_t child_cell_index1 = _grid.insert_cell_(cell_above_faces, tmp_faces, cell.marker());
  const size_t child_cell_index2 = _grid.insert_cell_(cell_below_faces, tmp_faces, cell.marker());
  _grid._n_inactive_cells++;
  if (_verbose)
  std::cout << "split " << parent_cell_index << " (parent "
            << _grid.cell(parent_cell_index).parent().index()
            << " ult " << _grid.cell(parent_cell_index).ultimate_parent().index() << "): "
            << child_cell_index1 << " " << child_cell_index2
            << std::endl << std::flush;


  // handle parent/child cell dependencies
  _grid.cell(parent_cell_index).m_children = {child_cell_index1, child_cell_index2};
  _grid.cell(child_cell_index1).m_parent = cell.index();
  _grid.cell(child_cell_index2).m_parent = cell.index();

  // we need to insert hanging nodes into neighboring cells
  // the eiasiest way to do it is to track split edges
  for (const auto & it_edge : find_affected_edges_(new_vertices, cell))
    for ( const auto icell : _grid.neighbors_indices_(it_edge.first) )
      if (icell != parent_cell_index && icell != child_cell_index1 && icell != child_cell_index2)
      {
        // std::cout << "insert hanging into " << _grid.cell(icell).index()
        //           << "(" << _grid.cell(icell).ultimate_parent().index() << ")"
        //           << std::endl;
        const auto & neighbor = _grid.cell(icell);
        if (!neighbor.has_vertex(it_edge.second))
          insert_hanging_node_(neighbor, it_edge.first, it_edge.second);
      }

  // finally, we need to split faces of the neighbors by face
  const auto & splitting_face = tmp_faces[split_face_local_index];
  for (std::size_t i=0; i<splitting_face.vertices.size(); ++i)
  {
    vertex_pair split_edge = {splitting_face.vertices[i], constants::invalid_index};
    if ( i + 1 <  splitting_face.vertices.size())
      split_edge.second = splitting_face.vertices[i+1];
    else split_edge.second = splitting_face.vertices[0];
    for ( const auto icell : _grid.neighbors_indices_(split_edge) )
      if (icell != cell.index() && icell != child_cell_index1 && icell != child_cell_index2)
      {
        const auto & neighbor = _grid.cell(icell);
        if (!neighbor.has_edge(split_edge))
        {
          // std::cout << "split neighbor face " << neighbor.index() << "("
          //           << neighbor.ultimate_parent().index() << std::endl;
          split_face_in_cell_(neighbor, split_edge);
        }

      }
  }
}

std::vector<std::size_t>
CellSplitter::build_global_face_indices_(const std::vector<size_t> & polygon_local_indices,
                                         const std::vector<size_t> & local_to_global) const
{
  std::vector<std::size_t> global_face_indices(polygon_local_indices.size());
  for (std::size_t i=0; i<polygon_local_indices.size(); ++i)
  {
    assert( polygon_local_indices[i] < local_to_global.size() );
    global_face_indices[i] = local_to_global[polygon_local_indices[i]];
  }

  return global_face_indices;
}

std::map<vertex_pair,size_t> CellSplitter::find_affected_edges_(const std::vector<size_t> &new_vertices,
                                                                const Cell & cell) const
{
  // map edge -> hanging vertex
  std::map<vertex_pair, size_t>  affected_edges;
  for (auto face : cell.faces())
    for (const vertex_pair & edge : face->edges())
      if ( affected_edges.find( std::minmax(edge.first, edge.second) ) == affected_edges.end() )
      {
        const Point & p1 = _grid.vertex(edge.first);
        const Point & p2 = _grid.vertex(edge.second);
        angem::Line<3, double> line(p1, p2-p1);
        for (const size_t v : new_vertices)
          if (line.distance(_grid.vertex(v)) < 1e-6)
            affected_edges.insert( {std::minmax(edge.first, edge.second), v} );
    }
  return affected_edges;
}

void CellSplitter::insert_hanging_node_(const Cell parent, const vertex_pair edge,
                                        const size_t inserted_vertex)
{
  assert( edge.first < edge.second );
  if (inserted_vertex == edge.first || inserted_vertex == edge.second)
    return;

  std::vector<FaceTmpData> tmp_faces;
  tmp_faces.reserve(parent.m_faces.size() + 2);

  for ( const auto face : parent.faces() )
  {
    // check whether the face contains the edge
    std::vector<vertex_pair> face_edges = face->edges();
    const size_t ne = face_edges.size();

    FaceTmpData f;
    f.parent = face->index();
    f.marker = face->marker();
    for (std::size_t i=0; i < ne; ++i)
    {
      if (i == 0)
        f.vertices.push_back(face_edges[i].first);
      if (std::min(face_edges[i].first,face_edges[i].second) == edge.first &&
          std::max(face_edges[i].first,face_edges[i].second) == edge.second)
       f.vertices.push_back(inserted_vertex);
      if (i != ne - 1)
        f.vertices.push_back(face_edges[i].second);
    }
    f.vtk_id = _grid.face_vtk_id_(f.vertices.size());
    tmp_faces.push_back( std::move(f) );
  }
  std::vector<size_t> indices_in_tmp(tmp_faces.size());
  std::iota(indices_in_tmp.begin(), indices_in_tmp.end(), 0);
  const size_t parent_cell_index = parent.index();
  const size_t child_cell_index = _grid.insert_cell_(indices_in_tmp, tmp_faces, parent.marker());
  _grid.cell(child_cell_index).m_parent = parent_cell_index;
  _grid.cell(parent_cell_index).m_children = {child_cell_index};
  _grid._n_inactive_cells++;
  if (_verbose)
  std::cout << "insert hanging into " << parent_cell_index << " "
            << "(" << _grid.cell(parent_cell_index).ultimate_parent().index() << "): "
            << child_cell_index
            << std::endl;
}

void CellSplitter::split_face_in_cell_(const Cell parent, const vertex_pair new_edge)
{
  std::vector<FaceTmpData> tmp_faces;
  tmp_faces.reserve(parent.m_faces.size() + 2);

  FaceTmpData f2;
  for (auto face : parent.faces())
  {
    FaceTmpData f1;
    f1.parent = face->index();
    f1.marker = face->marker();
    if (face->has_vertex(new_edge.first) && face->has_vertex(new_edge.second))
    {
      f2.parent = face->index();
      f2.marker = face->marker();
      bool in_first_part = true;
      for (const size_t v : face->vertices())
      {
        if (in_first_part) f1.vertices.push_back(v);
        else               f2.vertices.push_back(v);
        if (v == new_edge.first || v == new_edge.second)
        {
          in_first_part = !in_first_part;
          if (in_first_part) f1.vertices.push_back(v);
          else               f2.vertices.push_back(v);
        }
      }
      f1.vtk_id = _grid.face_vtk_id_(f1.vertices.size());
      f2.vtk_id = _grid.face_vtk_id_(f2.vertices.size());
      assert( f1.vertices.size() > 2 );
      assert( f2.vertices.size() > 2 );
    }
    else // leave as it is
    {
      f1.vertices = face->vertices();
      f1.vtk_id = face->vtk_id();
    }

    tmp_faces.push_back(std::move(f1));
  }
  if (!f2.vertices.empty())
  {
    tmp_faces.push_back(std::move(f2));
    std::vector<size_t> indices_in_tmp(tmp_faces.size());
    std::iota(indices_in_tmp.begin(), indices_in_tmp.end(), 0);
    const size_t parent_index = parent.index();
    const size_t child_cell_index = _grid.insert_cell_(indices_in_tmp, tmp_faces, parent.marker());
    _grid.cell(child_cell_index).m_parent = parent_index;
    _grid.cell(parent_index).m_children = {child_cell_index};
    _grid._n_inactive_cells++;
    if (_verbose)
    std::cout << "split neighbor face in "
              << parent_index << "(par "
              << _grid.cell(parent_index).parent().index() << ", ult "
              << _grid.cell(parent_index).ultimate_parent().index() << "): "
              << child_cell_index
              << std::endl;
  }
}


void CellSplitter::track_new_vertices_(std::vector<size_t> & global_vertex_indices,
                                       angem::PolyGroup<double> & split,
                                       std::vector<size_t> & new_vertices)
{
  for (size_t i = global_vertex_indices.size(); i < split.vertices.size(); ++i)
  {
    size_t new_vertex_index = _grid.n_vertices();
    // first check that the vertex isn't already present because it's been split
    const size_t split_v_index = _new_vertex_coord.find(split.vertices[i]);
    if (split_v_index == _new_vertices.size())
    {
      _new_vertex_coord.insert( split.vertices[i] );
      _new_vertices.push_back(new_vertex_index);
      new_vertices.push_back(new_vertex_index);
      _grid.insert_vertex(split.vertices[i]);
    }
    else
      new_vertex_index = _new_vertices[split_v_index];
    global_vertex_indices.push_back(new_vertex_index);
  }
}

void CellSplitter::create_face_groups_(angem::PolyGroup<double> & split,
                                       const std::vector<std::vector<size_t>> & face_vertex_global_numbering,
                                       const std::vector<size_t> & polygroup_polygon_parents,
                                       const std::vector<Face*> & cell_faces,
                                       const int splitting_face_marker,
                                       const size_t split_face_local_index,
                                       std::vector<size_t> & cell_above_faces,
                                       std::vector<size_t> & cell_below_faces,
                                       std::vector<FaceTmpData> & tmp_faces)
{
  tmp_faces.resize(split.polygons.size());
  for (size_t i = 0; i < split.polygons.size(); i++)
  {
    FaceTmpData & f = tmp_faces[i];
    f.vertices = face_vertex_global_numbering[i];
    f.vtk_id = _grid.face_vtk_id_(f.vertices.size());
    // determine the global index of the parent of the split face
    // splitting face is marked with polygroup_polygon_parents.size()
    f.parent = (polygroup_polygon_parents[i] < polygroup_polygon_parents.size()) ?
               cell_faces[ polygroup_polygon_parents[i] ]->index() :
               constants::invalid_index;

    if ( i == split_face_local_index )
      f.marker = splitting_face_marker;
    else
      f.marker = _grid.face(f.parent).marker();

    if ( split.markers[i] == constants::marker_below_splitting_plane ||
         split.markers[i] == constants::marker_splitting_plane )
      cell_below_faces.push_back(i);

    if (split.markers[i] == constants::marker_above_splitting_plane ||
        split.markers[i] == constants::marker_splitting_plane)
      cell_above_faces.push_back(i);
  }
}

}  // end namespace mesh
