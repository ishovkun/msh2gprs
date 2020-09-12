#include "CartesianMeshBuilder.hpp"

namespace mesh {

CartesianMeshBuilder::CartesianMeshBuilder(const CartesianMeshParameters & data)
    : _data(data)
{
  if (n_cells() == 0)
    throw std::invalid_argument("Cannot create grid with 0 cells");
  validate_stepping_(_data.dx);
  validate_stepping_(_data.dy);
  validate_stepping_(_data.dz);
}

CartesianMeshBuilder::operator Mesh() const
{
  Mesh grid;
  build(grid);
  return grid;
}

void CartesianMeshBuilder::build(Mesh & grid) const
{
  setup_vertices_(grid);
  setup_cells_(grid);
}

void CartesianMeshBuilder::validate_stepping_(const std::vector<double> &h) const
{
  if (std::any_of(h.begin(), h.end(), [](const double dx) {return std::fabs(dx) < 1e-10;}))
    throw std::invalid_argument("wrong stepping in cartesian grid");
}

size_t CartesianMeshBuilder::nx() const noexcept
{
  return _data.dx.size();
}

size_t CartesianMeshBuilder::ny() const noexcept
{
  return _data.dy.size();
}

size_t CartesianMeshBuilder::nz() const noexcept
{
  return _data.dz.size();
}

size_t CartesianMeshBuilder::nvx() const noexcept
{
  return nx() + 1;
}

size_t CartesianMeshBuilder::nvy() const noexcept
{
  return ny() + 1;
}

size_t CartesianMeshBuilder::nvz() const noexcept
{
  return nz() + 1;
}

size_t CartesianMeshBuilder::cell_index(size_t i, size_t j, size_t k) const
{
  if (i >= nx() || j >= ny() || k >= nz())
    throw std::invalid_argument("Wrong cell index");

  return ny()*nx()*k + nx()*j + i;
}

size_t CartesianMeshBuilder::vertex_index(size_t i, size_t j, size_t k) const
{
  if (k >= nvx() || j >= nvy() || k >= nvz())
    throw std::invalid_argument("Wrong vertex index");
  return nvy()*nvx()*k + nvx()*j + i;
}

size_t CartesianMeshBuilder::n_cells() const noexcept
{
  return nx()*ny()*nz();
}

size_t CartesianMeshBuilder::n_vertices() const noexcept
{
  return nvx() * nvy() * nvz();
}

void CartesianMeshBuilder::setup_vertices_(Mesh & grid) const
{
  grid.reserve_vertices(n_vertices());
  auto p = _data.origin;
  for (size_t k = 0; k < nvz(); ++k)
  {
    p[1] = _data.origin[1];
    for (size_t j = 0; j < nvy(); ++j)
    {
      p[0] = _data.origin[0];
      for (size_t i = 0; i < nvx(); ++i)
      {
        grid.insert_vertex(p);
        if (i < nx()) p[0] += _data.dx[i];
      }
      if (j < ny()) p[1] += _data.dy[j];
    }
    if (k < nz()) p[2] += _data.dz[k];
  }
  assert( grid.n_vertices() == n_vertices() );
}

void CartesianMeshBuilder::setup_cells_(Mesh & grid) const
{
  const int domain_marker = 0;
  const int left_boundary = 0;
  const int right_boundary = 1;
  const int front_boundary = 2;
  const int back_boundary = 3;
  const int bottom_boundary = 4;
  const int top_boundary = 5;

  grid.cells().reserve(n_cells());
  for (size_t k = 0; k < nz(); ++k)
    for (size_t j = 0; j < ny(); ++j)
      for (size_t i = 0; i < nx(); ++i)
      {
        // create a hex: VTK vertex numbering
        const size_t v0 = vertex_index(i,   j,   k);
        const size_t v1 = vertex_index(i+1, j,   k);
        const size_t v2 = vertex_index(i+1, j+1, k);
        const size_t v3 = vertex_index(i,   j+1, k);
        const size_t v4 = vertex_index(i,   j,   k+1);
        const size_t v5 = vertex_index(i+1, j,   k+1);
        const size_t v6 = vertex_index(i+1, j+1, k+1);
        const size_t v7 = vertex_index(i,   j+1, k+1);
        const std::vector<size_t> cell_vertices = {v0, v1, v2, v3, v4, v5, v6, v7};
        grid.insert_cell(cell_vertices, angem::HexahedronID, domain_marker );
        if (i == 0)
          grid.insert_face({v0, v4, v7, v3}, angem::QuadrangleID, left_boundary);
        else if (i == nx() - 1)
          grid.insert_face({v1, v2, v6, v5}, angem::QuadrangleID, right_boundary);
        if (j == 0)
          grid.insert_face({v0, v1, v5, v4}, angem::QuadrangleID, front_boundary);
        else if (j == ny() - 1)
          grid.insert_face({v2, v3, v7, v6}, angem::QuadrangleID, back_boundary);
        if (k == 0)
          grid.insert_face({v0, v1, v2, v3}, angem::QuadrangleID, bottom_boundary);
        else if (k == nz() - 1)
          grid.insert_face({v4, v7, v6, v5}, angem::QuadrangleID, top_boundary);
      }
}

}  // end namespace mesh
