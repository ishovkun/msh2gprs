#include "PolyhedralElementBase.hpp"
#include "gmsh_interface/GmshInterface.hpp"  // provides gprs_data::GmshInterface
#include "mesh/Subdivision.hpp"              // provides mesh::Subdivision
#include "VTKWriter.hpp"                     // provides VTKWriter

#ifdef WITH_EIGEN

namespace discretization {

using api = gprs_data::GmshInterface;

PolyhedralElementBase::PolyhedralElementBase(const mesh::Cell & cell,
                                             const mesh::Mesh & grid,
                                             const FiniteElementConfig & config)
    : _parent_cell(cell), _parent_grid(grid), _config(config)
{}

void PolyhedralElementBase::build_triangulation_()
{
  // triangulate the polyhedral element
  if (_config.subdivision_method == PolyhedralFEMSubdivision::gmsh_generate)
  {
    api::initialize_gmsh();
    api::build_triangulation(_parent_cell, _element_grid, double(_config.order));
    api::finalize_gmsh();
  }
  else if (_config.subdivision_method == PolyhedralFEMSubdivision::refinement)
  {
    mesh::Subdivision subdivision(_parent_cell, _element_grid, _config.order);
  }
  else throw std::invalid_argument("unknown subdivision method");
}

std::vector<std::vector<size_t>> PolyhedralElementBase::create_face_domains_()
{
  std::vector<std::vector<size_t>> parent_face_children(_parent_cell.faces().size());
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const size_t parent_face_index = static_cast<size_t>(face->marker() - 1);
      parent_face_children[ parent_face_index ].push_back( face->index() );
    }
  return parent_face_children;
}

void PolyhedralElementBase::save_shape_functions(const std::string fname) const
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _parent_cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0.0);
    for (size_t i = 0; i < nv; ++i)
      output[i] = _basis_functions[j][i];
    IO::VTKWriter::add_data(output, "basis-"+std::to_string(j), out);
  }
  out.close();
}


}  // end namespace discretization

#endif
