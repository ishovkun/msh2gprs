#include "PolyhedralElementBase.hpp"
#include "gmsh_interface/GmshInterface.hpp"  // provides gprs_data::GmshInterface
#include "mesh/Subdivision.hpp"              // provides mesh::Subdivision
#include "mesh/io/VTKWriter.hpp"                     // provides VTKWriter
#include "PFEM_integration/TributaryRegion2dFaces.hpp"
#include "PFEM_integration/TributaryRegion3dFaces.hpp"
#include "PFEM_integration/TributaryRegion2dVertices.hpp"
#include "PFEM_integration/TributaryRegion3dVertices.hpp"
#include "PFEM_integration/TributaryRegion3dFull.hpp"
#include "PFEM_integration/IntegrationRule3d.hpp"
#include "PFEM_integration/IntegrationRule2d.hpp"
#include "PFEM_integration/IntegrationRuleFracture.hpp"
// #include "PFEM_integration/IntegrationRuleFractureAverage.hpp"  // provides IntegrationFractureAverage
// #include "PFEM_integration/IntegrationRuleFractureFull.hpp"  // provides IntegrationFractureFull
#include "PFEM_integration/FaceSorter.hpp"  // provides Facesorter
#include <memory>


#ifdef WITH_EIGEN

namespace discretization {

using api = gprs_data::GmshInterface;

PolyhedralElementBase::PolyhedralElementBase(const mesh::Cell & cell,
                                             const mesh::Mesh & grid,
                                             const FiniteElementConfig & config,
                                             const bool sort_faces)
    : _parent_cell(cell), _parent_grid(grid), _config(config), _sort_faces(sort_faces)
{}

void PolyhedralElementBase::build_triangulation_()
{
  // triangulate the polyhedral element
  if (_config.subdivision_method == PolyhedralFEMSubdivision::gmsh_generate)
  {
    api::initialize_gmsh();
    api::build_triangulation(_parent_cell, _subgrid, double(_config.order));
    api::finalize_gmsh();
  }
  else if (_config.subdivision_method == PolyhedralFEMSubdivision::refinement)
  {
    mesh::Subdivision subdivision(_parent_cell, _subgrid, _config.order);
  }
  else throw std::invalid_argument("unknown subdivision method");
}

std::vector<std::vector<size_t>> const & PolyhedralElementBase::get_face_domains()
{
  if (_face_domains.empty())
    _face_domains = create_face_domains_();
  return _face_domains;
}

std::vector<std::vector<size_t>> PolyhedralElementBase::create_face_domains_()
{
  std::vector<std::vector<size_t>> parent_face_children(_parent_cell.faces().size());
  for (auto face = _subgrid.begin_active_faces(); face != _subgrid.end_active_faces(); ++face)
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const size_t parent_face_index = static_cast<size_t>(face->marker() - 1);
      parent_face_children[ parent_face_index ].push_back( face->index() );
    }

  // the order of fe_face_values that correspond to fractures must be consistent with the
  // neighbor's fe_face_values
  if (_sort_faces)
    for (size_t iface = 0; iface < parent_face_children.size(); ++iface)
    {
      FaceSorter sorter(*_parent_cell.faces()[iface], _subgrid);
      sorter.sort(parent_face_children[iface]);
    }

  return parent_face_children;
}

void PolyhedralElementBase::save_shape_functions(const std::string fname) const
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  mesh::IO::VTKWriter::write_geometry(_subgrid, out);
  const size_t nv = _subgrid.n_vertices();
  mesh::IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _parent_cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0.0);
    for (size_t i = 0; i < nv; ++i)
      output[i] = _basis_functions[j][i];
    mesh::IO::VTKWriter::add_data(output, "basis-"+std::to_string(j), out);
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
  std::unique_ptr<TributaryRegion3dBase> regions = nullptr;
  switch (_config.integration_rule)
  {
    case PolyhedronIntegrationRule::Full:
      {
        regions = std::make_unique<TributaryRegion3dFull>(*this);
        break;
      }
    case PolyhedronIntegrationRule::FacesAverage:
      {
        // TributaryRegion3dFaces tributary3d(*this);
        regions = std::make_unique<TributaryRegion3dFaces>(*this);
        // IntegrationRule3dAverage rule_cell(*this, tributary3d);
        break;
      }
    // case PolyhedronIntegrationRule::FacesPointwise:
    //   {
    //     regions = std::make_unique<TributaryRegion3dFaces>(*this);
    //     TributaryRegion3dFaces tributary3d(*this);
    //     IntegrationRule3dPointwise rule_cell(*this, tributary3d);
    //     break;
    //   }
    case PolyhedronIntegrationRule::VerticesAverage:
      {
        regions = std::make_unique<TributaryRegion3dVertices>(*this);
        break;
      }
    // case PolyhedronIntegrationRule::VerticesPointwise:
    //   {
    //     TributaryRegion3dVertices tributary3d(*this);
    //     IntegrationRule3dAverage rule_cell(*this, tributary3d);
    //     break;
    //   }
    default:
      throw std::invalid_argument("Integration rule unknown");
  }

  _integration_rule3d = std::make_shared<IntegrationRule3d>(*this, *regions);
  _cell_data = _integration_rule3d->integrate(_parent_cell.vertex_coordinates());
}

FiniteElementData PolyhedralElementBase::get_face_data(size_t iface)
{
  auto const faces = _parent_cell.faces();
  auto const & rule = integration_rule2(iface);
  auto const basis = get_face_basis_(*_parent_cell.faces()[iface], _parent_cell);
  return rule.integrate(faces[iface]->vertex_coordinates(), basis);
}

FiniteElementData PolyhedralElementBase::get_fracture_data(const size_t iface,
                                                           const angem::Basis<3,double> basis)
{
  // auto const basis1 = get_face_basis_(*_parent_cell.faces()[iface], _parent_cell);
  return integration_rule_frac(iface).integrate(_parent_cell.vertex_coordinates(), basis);
}

std::unique_ptr<TributaryRegion2dBase> PolyhedralElementBase::build_tributary_2d_(const size_t parent_face)
{
  switch (_config.integration_rule)
  {
    // case PolyhedronIntegrationRule::VerticesPointwise:
    //   {
    //     TributaryRegion2dVertices tributary2d(*this);
    //     _tributary_2d = tributary2d.get();
    //   break;
    // }
    case PolyhedronIntegrationRule::FacesAverage:
      return std::make_unique<TributaryRegion2dFaces>(*this, parent_face);
    case PolyhedronIntegrationRule::VerticesAverage:
      return std::make_unique<TributaryRegion2dVertices>(*this, parent_face);
    // case PolyhedronIntegrationRule::Full:
    //   return nullptr;
    //   break;
    default:
      throw std::invalid_argument("Tributary region not implemented");
  }
}

IntegrationRule2d const & PolyhedralElementBase::integration_rule2(size_t iface)
{
  if ( _integration_rules2d.empty() )  // lazy initialize
    _integration_rules2d.resize(_parent_cell.faces().size(), nullptr);

  if (!_integration_rules2d[iface])  // lazy evaluate and store
  {
    auto const tributary_region = build_tributary_2d_(iface);
    auto const basis = get_face_basis_(*_parent_cell.faces()[iface], _parent_cell);
    _integration_rules2d[iface] = std::make_shared<IntegrationRule2d>(*this, *tributary_region, iface, basis);
  }
  return *_integration_rules2d[iface];
}

IntegrationRuleFracture const & PolyhedralElementBase::integration_rule_frac(size_t iface)
{
  if ( _integration_rules_frac.empty() )  // lazy initialize
    _integration_rules_frac.resize(_parent_cell.faces().size(), nullptr);

  if (  !_integration_rules_frac[iface]  )
  {
    auto const trib = build_tributary_2d_(iface);
    auto const basis = get_face_basis_(*_parent_cell.faces()[iface], _parent_cell);
    _integration_rules_frac[iface] = std::make_shared<IntegrationRuleFracture>(*this, *trib, iface, basis);
  }
  return *_integration_rules_frac[iface];
}


}  // end namespace discretization

#endif
