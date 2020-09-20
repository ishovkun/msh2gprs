#include "RefinementAspectRatio.hpp"
#include "io/VTKWriter.hpp"

namespace mesh {

RefinementAspectRatio::RefinementAspectRatio(Mesh & grid, const double aspect_ratio, const size_t max_level)
    : _grid(grid), _ratio(aspect_ratio), _max_level(max_level), _splitter(_grid)
{
  if (aspect_ratio < 1)
    throw std::invalid_argument("aspect ratio must be > 1");

  split_cells_();
}

void RefinementAspectRatio::split_cells_()
{
  std::list<size_t> flag_refine;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    if (check_cell_(cell->index()))
      flag_refine.push_back(cell->index());

  for (const size_t cell_id : flag_refine)
  {
    std::cout << "cell_id = " << cell_id << std::endl;
    const auto cut_plane = find_cut_plane_(cell_id);
    std::vector<angem::Point<3,double>> points = {
     cut_plane.origin(),
     cut_plane.origin() + 10 * cut_plane.get_basis()(0) ,
     cut_plane.origin() + 10 * cut_plane.get_basis()(1) ,
    };
    // angem::Polygon<double> poly (points);
    // IO::VTKWriter::write<double>(poly, "mycut.vtk");
    // IO::VTKWriter::write_geometry(_grid, _grid.cell(cell_id),
                                  // "cell-"+std::to_string(cell_id) + ".vtk");

    const auto & cell = _grid.cell(cell_id);
    if (cell.level() < _max_level)
      _splitter.split_cell(_grid.cell(cell_id), cut_plane, constants::default_face_marker);
    // IO::VTKWriter::write_geometry(_grid, _grid.cell(20), "cell-20.vtk");
    // IO::VTKWriter::write_geometry(_grid, _grid.cell(35), "cell-35.vtk");
    // IO::VTKWriter::write_geometry(_grid, _grid.cell(36), "cell-36.vtk");
    // exit(0);
    // return;
  }
}

bool RefinementAspectRatio::check_cell_(const size_t cell_id) const
{
  const auto & cell = _grid.cell(cell_id);
  double min_len = std::numeric_limits<double>::max();
  double max_len = std::numeric_limits<double>::lowest();
  for (const auto & edge : cell.edges())
  {
    const double l = _grid.vertex(edge.first).distance(_grid.vertex(edge.second));
    min_len = std::min(l, min_len);
    max_len = std::max(l, max_len);
  }

  return  max_len / min_len > _ratio;
}

angem::Plane<double> RefinementAspectRatio::find_cut_plane_(const size_t cell_id) const
{
  const auto & cell = _grid.cell(cell_id);
  auto edges = cell.edges();
  std::sort(edges.begin(), edges.end(),
            [this](const auto & edge1, const auto & edge2) {
              const double l1 = _grid.vertex(edge1.first).distance(_grid.vertex(edge1.second));
              const double l2 = _grid.vertex(edge2.first).distance(_grid.vertex(edge2.second));
              return l1 < l2;
            });
  const auto & e1 = edges.back();
  const auto & e2 = edges[edges.size() - 2];
  const auto c1 = edge_center_(e1);
  const auto c2 = edge_center_(e2);

  const auto t1 = c2 - c1;
  const auto t2 = _grid.vertex(e1.second) - _grid.vertex(e1.first);
  angem::Plane<double> pln1(c1, c2, _grid.vertex(e1.second));
  angem::Point<3,double> c3;
  bool found = false;
  for (int i = edges.size() - 3; i > 0; i--)
  {
    c3 = edge_center_(edges[i]);
    if (std::fabs(pln1.signed_distance(c3)) > 1e-3*edge_length_(edges[i]))
    {
      found = true;
      break;
    }
  }
  if (found) return angem::Plane<double>(c1, c2, c3);
  else throw std::runtime_error("could not find a feasible cut to split");
}

angem::Point<3,double> RefinementAspectRatio::edge_center_(const vertex_pair & edge) const
{
  return 0.5*(_grid.vertex(edge.first) + _grid.vertex(edge.second));
}

double RefinementAspectRatio::edge_length_(const vertex_pair & edge) const
{
  return _grid.vertex(edge.first).distance(_grid.vertex(edge.second));
}
}  // end namespace mesh
