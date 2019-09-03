// #include <face_iterator.hpp>
// #include <mesh_methods.hpp>

// namespace mesh
// {

// face_iterator::
// face_iterator(const std::size_t                       face_index,
//               std::vector<Face>                     & grid_faces,
//               std::vector<angem::Point<3,double>>   & grid_vertices)
//     : face_index(face_index),
//       grid_faces(grid_faces),
//       grid_vertices(grid_vertices)
// {}


// face_iterator::
// face_iterator(const face_iterator & other)
//     : face_index(other.face_index),
//       grid_faces(other.grid_faces),
//       grid_vertices(other.grid_vertices)
// {}


// bool face_iterator::operator==(const face_iterator & other) const
// {
//   if (face_index != other.face_index)
//     return false;
//   else
//     return true;
// }


// bool face_iterator::operator!=(const face_iterator & other) const
// {
//   return !(*this == other);
// }


// face_iterator & face_iterator::operator++()
// {
//   face_index++;
//   if (face_index > grid_faces.size())
//     std::out_of_range("face iterator out of range");
//   return (*this);
// }


// int face_iterator::marker()
// {
//   return operator*()->marker;
// }


// std::vector<Point> face_iterator::vertex_coordinates() const
// {
//   return get_vertex_coordinates(grid_vertices, operator*()->vertices);
// }


// std::vector<std::size_t> face_iterator::vertex_indices() const
// {
//   return operator*()->vertices;
// }


// Point face_iterator::normal() const
// {
//   const auto poly = angem::Polygon<double>(vertex_coordinates());
//   return poly.plane.normal();
// }


// angem::Polygon<double> face_iterator::polygon() const
// {
//   return angem::Polygon<double>(vertex_coordinates());
// }


// angem::Point<3,double> face_iterator::center() const
// {
//   const auto verts = vertex_indices();
//   angem::Point<3,double> c = {0, 0, 0};
//   for (const size_t vert : verts)
//     c += grid_vertices[vert];

//   for (int comp = 0; comp < 3; comp++)
//     c[comp] /= verts.size();
//   return c;
// }

// }
