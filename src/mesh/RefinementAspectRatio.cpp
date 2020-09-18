#include "RefinementAspectRatio.hpp"

namespace mesh {

RefinementAspectRatio::RefinementAspectRatio(Mesh & grid, const double aspect_ratio, const size_t max_level)
    : _grid(grid), _ratio(aspect_ratio), _max_level(max_level)
{
  if (aspect_ratio < 1)
    throw std::invalid_argument("aspect ratio must be > 1");

  find_problematic_cells_();
}

void RefinementAspectRatio::find_problematic_cells_()
{
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    check_cell_(cell->index());
  exit(0);
}

bool RefinementAspectRatio::check_cell_(const size_t cell_id) const
{
  const auto & cell = _grid.cell(cell_id);
  const auto edges = cell.edges();
  double min_len = std::numeric_limits<double>::max();
  double max_len = std::numeric_limits<double>::lowest();
  for (const auto & edge : cell.edges())
  {
    const double l = _grid.vertex(edge.first).distance(_grid.vertex(edge.second));
    min_len = std::min(l, min_len);
    max_len = std::max(l, max_len);
  }

  std::cout << cell_id << ": ratio = " << max_len / min_len << std::endl;
  if ( max_len / min_len > _ratio )
  {
    std::cout << "found" << std::endl;
    std::cout << "cell_id = " << cell_id << std::endl;
    exit(0);
    return true;
  }
  else return false;
}

}  // end namespace mesh
