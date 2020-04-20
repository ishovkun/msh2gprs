#include "PolyhedralElementBase.hpp"
#include "gmsh_interface/GmshInterface.hpp"  // provides gprs_data::GmshInterface
#include "mesh/Subdivision.hpp"              // provides mesh::Subdivision

namespace discretization {

using api = gprs_data::GmshInterface;

PolyhedralElementBase::PolyhedralElementBase(const mesh::Cell & cell,
                                             const FiniteElementConfig & config)
    : _parent_cell(cell), _config(config)
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
  //   std::string fname = "custom_subdivision.vtk";
  //   std::cout << "saving " << fname << std::endl;
  //   std::ofstream out;
  //   out.open(fname.c_str());
  //   IO::VTKWriter::write_geometry(_element_grid, out);
  // IO::VTKWriter::enter_section_point_data(_element_grid.n_vertices(), out);
  // std::vector<double> output(_element_grid.n_vertices(), 0);
  // for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
  // {
  //   if (face->marker() > 0 && face->neighbors().size() == 1)
  //   {
  //     for ( const size_t v : face->vertices() )
  //       output[v] = face->marker();
  //   }
  // }
  // IO::VTKWriter::add_data(output, "bnd-marker", out);
    // exit(0);
  }
  else throw std::invalid_argument("unknown subdivision method");
}

std::vector<std::vector<size_t>> PolyhedralElementBase::create_face_domains_()
{
  std::vector<std::vector<size_t>> parent_face_children(_parent_cell.faces().size());
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const size_t parent_face_index = face->marker() - 1;
      parent_face_children[ parent_face_index ].push_back( face->index() );
    }

  return parent_face_children;
}

}  // end namespace discretization
