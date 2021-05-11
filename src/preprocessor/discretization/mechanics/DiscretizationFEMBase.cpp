#include "DiscretizationFEMBase.hpp"
#include "logger/ProgressBar.hpp"  // provides ProgressBar
#include "GlobalOpts.hpp"

namespace discretization {

DiscretizationFEMBase::DiscretizationFEMBase(mesh::Mesh & grid, const FiniteElementConfig & config,
                                             const std::vector<int> & fracture_face_markers,
                                             const std::vector<size_t> & neumann)
    : _grid(grid), _config( config ),
      _neumann_faces(neumann.begin(), neumann.end())
{
  for (const int marker : fracture_face_markers)
    _fracture_face_orientation[marker] = {angem::Basis<3, double>(), false};

  _face_data.resize( _grid.n_faces() );
  _cell_data.resize( _grid.n_cells_total() );
  _frac_data.resize( _grid.n_faces() );
}

void DiscretizationFEMBase::build()
{
  std::unique_ptr<logging::ProgressBar> progress = nullptr;
  if (gprs_data::GlobalOpts::ref().print_progressbar)
    progress = std::make_unique<logging::ProgressBar>("Build Finite Elements", _grid.n_active_cells());

  for (auto [cell, counter] = std::tuple( _grid.begin_active_cells(), 0 ); cell != _grid.end_active_cells(); ++cell, ++counter)
  {
    if (progress) progress->set_progress(counter);
    // std::cout << "building cell " << cell->index() << std::endl;
    build_(*cell);
    build_cell_data_(*cell);
    build_face_data_(*cell);
  }

  if (progress) progress->finalize();
}

void DiscretizationFEMBase::build_cell_data_(mesh::Cell const & cell)
{
  FiniteElementData cell_fem_data = _element->get_cell_data();
  cell_fem_data.element_index = cell.index();
  _cell_data[cell.index()] = std::move(cell_fem_data);
}

void DiscretizationFEMBase::build_face_data_(mesh::Cell const & cell)
{
  auto const faces = cell.faces();
  for ( size_t f = 0; f < faces.size(); ++f)
  {
    auto const * face = faces[f];
    const size_t face_index = face->index();
    const bool is_fracture = _fracture_face_orientation.find( face->marker() ) != _fracture_face_orientation.end();
    const bool is_neumann = _neumann_faces.find( face->index() ) != _neumann_faces.end();
    if (is_fracture)
      _face_basis = get_basis_(*face, _fracture_face_orientation.find(face->marker())->second);

    if ( (is_neumann || is_fracture) && _face_data[face->index()].points.empty() )
    {
      _face_data[face_index] = _element->get_face_data(f);
      _face_data[face_index].element_index = face_index;
    }

    if (is_fracture)
    {
      _frac_data[face_index].push_back(_element->get_fracture_data(f, _face_basis));
      _frac_data[face_index].back().element_index = cell.index();
    }
  }
}

const angem::Basis<3, double> &
DiscretizationFEMBase::get_basis_(const mesh::Face & face,
                                  FaceOrientation &orientation) noexcept
{
  if (!orientation.assigned)
  {
    orientation.basis = face.polygon().plane().get_basis();
    orientation.assigned = true;
  }

  return orientation.basis;
}


}  // end namespace discretization
