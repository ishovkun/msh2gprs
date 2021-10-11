#include "DiscretizationINSIM.hpp"

namespace discretization {

DiscretizationINSIM::DiscretizationINSIM(DoFNumbering const & dof_numbering,
                                         gprs_data::SimData & data,
                                         std::vector<ControlVolumeData> & cv_data,
                                         std::vector<ConnectionData> & connection_data)
    : DiscretizationBase(dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationINSIM::build()
{
  for (size_t v = 0; v < m_grid.n_vertices(); ++v)
    if ( m_dofs.has_vertex( v ) )
      build_vertex_data_( v );

}

void DiscretizationINSIM::build_vertex_data_(size_t vertex)
{
  auto const attached_cells = m_grid.vertex_cells( vertex );
  size_t const n = attached_cells.size();

  std::vector<ControlVolumeData> attached_cvs( n );
  for (size_t i = 0; i < n; ++i)
    DiscretizationBase::build_cell_data_( *attached_cells[i], attached_cvs[i] );

  size_t const dof = m_dofs.vertex_dof( vertex );

  auto & cv = m_cv_data[dof];
  cv.type = ControlVolumeType::cell;
  cv.master = vertex;
  cv.center = m_grid.vertex( vertex );

  // take average of physical properties
  cv.custom.resize( attached_cvs.front().custom.size(), 0.f );
  for (size_t i = 0; i < n; ++i) {
    cv.volume += attached_cvs[i].volume;
    cv.porosity += attached_cvs[i].porosity / n;
    cv.permeability += (attached_cvs[i].permeability / (double)n);
    for (size_t j = 0; j < cv.custom.size(); ++j)
      cv.custom[j] += attached_cvs[i].custom[j] / n;
  }
}


}  // end namespace discretization
