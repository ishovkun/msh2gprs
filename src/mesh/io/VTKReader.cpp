#include "VTKReader.hpp"

namespace mesh {

namespace io {

VTKReader::VTKReader(const std::string file_name, mesh::Mesh & grid)
    : _grid(grid)
{
  read_file_(file_name);
}

void VTKReader::read_file_(const std::string & file_name)
{
  std::fstream in;
  in.open(file_name.c_str(), std::fstream::in);
  if (!in) throw std::out_of_range(file_name + " does not exist");

  read_header_(in);
  read_vertices_(in);
  read_cells_(in);
  read_cell_types_(in);
  in.close();
  create_grid_();
}

void VTKReader::read_header_(std::fstream & in) const
{
  std::string line;
  std::getline(in, line);  // vtk DataFile Version 3.0
  std::getline(in, line);  // 3D Grid
  std::getline(in, line);  // ASCII
  if (line != "ASCII")
    throw std::invalid_argument("File format not supported");
  std::getline(in, line);  // DATASET UNSTRUCTURED_GRID
}

void VTKReader::read_vertices_(std::fstream & in)
{
  std::string s;
  in >> s;
  if (s != "POINTS")
    throw std::invalid_argument("Vertices should start with POINTS");

  size_t n_vertices;
  in >> n_vertices;
  if (n_vertices == 0)
    throw std::invalid_argument("POINTS should have more than 0 vertices");
  std::cout << "n_vertices = " << n_vertices << std::endl;

  in >> s;
  if (s != "float")
    throw std::invalid_argument("Unknown POINTS format (should be FLOAT)");

  auto & vertex_coord = _grid.vertices();
  vertex_coord.resize(n_vertices);
  for (size_t i=0; i<n_vertices; ++i)
    for (size_t j=0; j<3; ++j)
      in >> vertex_coord[i][j];
}

void VTKReader::read_cells_(std::fstream & in)
{
  std::string s;
  in >> s;
  if ( s != "CELLS" )
    throw std::invalid_argument("Entry should be CELLS");
  size_t n_cells;
  in >> n_cells;
  std::cout << "n_cells = " << n_cells << std::endl;

  size_t n_cell_entries;
  in >> n_cell_entries;
  _cell_entries.resize(n_cell_entries);

  for (size_t i=0; i<n_cell_entries; ++i)
    in >> _cell_entries[i];
}

void VTKReader::read_cell_types_(std::fstream & in)
{
  std::string s;
  in >> s;
  if (s != "CELL_TYPES")
  {
    std::cout << s << std::endl;
    throw std::invalid_argument("Entry should be CELL_TYPES");
  }

  size_t n_cell_type_entries;
  in >> n_cell_type_entries;
  _vtk_ids.resize(n_cell_type_entries);
  for (size_t i=0; i<n_cell_type_entries; ++i)
    in >> _vtk_ids[i];
}

void VTKReader::create_grid_()
{
  const size_t n_cells = _vtk_ids.size();
  size_t idx = 0;
  for (size_t icell=0; icell < n_cells; ++icell)
  {
    const int id = _vtk_ids[icell];
    if (id == angem::VTK_ID::GeneralPolyhedronID)
      create_general_polyhedron_cell_(idx);
    else
      create_regular_polyhedron_cell_(id, idx);
  }
}

void VTKReader::create_general_polyhedron_cell_(size_t & idx)
{
  const size_t n_cell_entries = _cell_entries[idx++];
  const size_t n_faces = _cell_entries[idx++];
  std::vector<mesh::FaceTmpData> faces(n_faces);
  for (size_t iface=0; iface<n_faces; ++iface)
  {
    const size_t nv = _cell_entries[idx++];
    auto & face = faces[iface];
    face.vertices.resize(nv);
    for (size_t i=0; i<nv; ++i)
      face.vertices[i] = _cell_entries[idx++];
  }

  std::vector<size_t> take_faces(n_faces);
  std::iota(take_faces.begin(), take_faces.end(), 0);
  std::set<size_t> cell_vertices_set;
  for (auto & f : faces)
    for (const size_t v : f.vertices)
      cell_vertices_set.insert(v);
  std::vector<size_t> cell_vertices(cell_vertices_set.begin(), cell_vertices_set.end());


  _grid.insert_cell_( cell_vertices, take_faces, faces , angem::VTK_ID::GeneralPolyhedronID,
                      mesh::constants::default_cell_marker);
}

void VTKReader::create_regular_polyhedron_cell_(const int id, size_t & idx)
{
  const size_t n_cell_entries = _cell_entries[idx++];
  std::vector<size_t> cell_vertices(n_cell_entries);
  for (size_t i=0; i<n_cell_entries; ++i)
    cell_vertices[i] = _cell_entries[idx++];
  _grid.insert_cell(cell_vertices, id);
}

}  // end namespace io

}  // end namespace mesh
