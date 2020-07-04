#include "PolyhedralElementBase.hpp"
#include "gmsh_interface/GmshInterface.hpp"  // provides gprs_data::GmshInterface
#include "mesh/Subdivision.hpp"              // provides mesh::Subdivision
#include "VTKWriter.hpp"                     // provides VTKWriter
#include "PFEM_integration/TributaryRegion2dFaces.hpp"
#include "PFEM_integration/TributaryRegion3dFaces.hpp"
#include "PFEM_integration/TributaryRegion2dVertices.hpp"
#include "PFEM_integration/TributaryRegion3dVertices.hpp"
#include "PFEM_integration/IntegrationRule3dAverage.hpp"
#include "PFEM_integration/IntegrationRule3dFull.hpp"
#include "PFEM_integration/IntegrationRule3dPointwise.hpp"
#include "PFEM_integration/IntegrationRule2dAverage.hpp"
#include "PFEM_integration/IntegrationRule2dPointwise.hpp"
#include "PFEM_integration/IntegrationRule2dFull.hpp"
#include "PFEM_integration/IntegrationRuleFractureAverage.hpp"  // provides IntegrationFractureAverage
#include "PFEM_integration/IntegrationRuleFractureFull.hpp"  // provides IntegrationFractureFull
#include "PFEM_integration/FaceSorter.hpp"  // provides Facesorter


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

std::vector<std::vector<size_t>> PolyhedralElementBase::create_face_domains_(const bool sort_faces)
{
  std::vector<std::vector<size_t>> parent_face_children(_parent_cell.faces().size());
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const size_t parent_face_index = static_cast<size_t>(face->marker() - 1);
      parent_face_children[ parent_face_index ].push_back( face->index() );
    }

  // the order of fe_face_values that correspond to fractures must be consistent with the
  // neighbor's fe_face_values
  if (sort_faces)
    for (size_t iface = 0; iface < parent_face_children.size(); ++iface)
    {
      FaceSorter sorter(*_parent_cell.faces()[iface], _element_grid);
      sorter.sort(parent_face_children[iface]);
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

std::vector<std::list<size_t>> PolyhedralElementBase::map_parent_vertices_to_parent_faces_()
{
  // map parent vertex to parent faces
  std::vector<std::list<size_t>> parent_vertex_markers( _parent_cell.n_vertices() );
  const std::vector<size_t> parent_vertices = _parent_cell.vertices();
  const std::vector<const mesh::Face*> parent_faces = _parent_cell.faces();
  for (size_t ipf=0; ipf<parent_faces.size(); ++ipf)
  {
    const auto parent_face = parent_faces[ipf];

    for (size_t iv_parent=0; iv_parent<parent_vertices.size(); ++iv_parent)
    {
      const size_t parent_vertex = parent_vertices[iv_parent];
      if (parent_face->has_vertex(parent_vertex))
        parent_vertex_markers[iv_parent].push_back(ipf + 1);
    }
  }
  return parent_vertex_markers;
}

void PolyhedralElementBase::build_fe_cell_data_()
{
  switch (_config.integration_rule)
  {
    case PolyhedronIntegrationRule::Full:
      {
        IntegrationRule3dFull rule_cell(*this);
        break;
      }
    case PolyhedronIntegrationRule::FacesAverage:
      {
        TributaryRegion3dFaces tributary3d(*this);
        IntegrationRule3dAverage rule_cell(*this, tributary3d);
        break;
      }
    case PolyhedronIntegrationRule::FacesPointwise:
      {
        TributaryRegion3dFaces tributary3d(*this);
        IntegrationRule3dPointwise rule_cell(*this, tributary3d);
        break;
      }
    case PolyhedronIntegrationRule::VerticesAverage:
      {
        TributaryRegion3dVertices tributary3d(*this);
        IntegrationRule3dAverage rule_cell(*this, tributary3d);
        break;
      }
    case PolyhedronIntegrationRule::VerticesPointwise:
      {
        TributaryRegion3dVertices tributary3d(*this);
        IntegrationRule3dAverage rule_cell(*this, tributary3d);
        break;
      }
    default:
      throw std::invalid_argument("Integration rule unknown");
  }
}

void PolyhedralElementBase::build_fe_face_data_()
{
  switch (_config.integration_rule)
  {
    case PolyhedronIntegrationRule::VerticesPointwise:
      {
        TributaryRegion2dVertices tributary2d(*this);
        IntegrationRule2dPointwise rule_faces(*this, tributary2d);
        break;
      }
    case PolyhedronIntegrationRule::FacesAverage:
      {
        TributaryRegion2dFaces tributary2d(*this);
        IntegrationRule2dAverage rule_faces(*this, tributary2d);
        break;
      }
    case PolyhedronIntegrationRule::Full:
      {
        IntegrationRule2dFull tributary2d(*this);
        break;
      }
    default:
      {
        TributaryRegion2dVertices tributary2d(*this);
        IntegrationRule2dPointwise rule_faces(*this, tributary2d);
        break;
      }
  }
}

void PolyhedralElementBase::build_fe_fracture_data_()
{
  switch (_config.integration_rule)
  {
    case PolyhedronIntegrationRule::Full:
      {
        IntegrationRuleFractureFull rule_fractures(*this);
        break;
      }
    case PolyhedronIntegrationRule::FacesAverage:
      {
        TributaryRegion2dFaces tributary2d(*this);
        IntegrationRuleFractureAverage rule_fractures(*this, tributary2d);
        break;
      }
    default:
      throw std::invalid_argument("Not implemented");
  }
}

}  // end namespace discretization

#endif
