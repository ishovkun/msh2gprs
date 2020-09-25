#include "GridIntersectionSearcher.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include "angem/Collisions.hpp"
#include "Logger.hpp"

namespace gprs_data {

GridIntersectionSearcher::GridIntersectionSearcher(const mesh::Mesh & grid)
    : _grid(grid), _mapper(grid)
{}

double GridIntersectionSearcher::top() const noexcept
{
  return get_highest_bound_(2);
}

double GridIntersectionSearcher::bottom() const noexcept
{
  return get_lowest_bound_(2);
}

double GridIntersectionSearcher::get_lowest_bound_(const size_t direction) const
{
  const auto & sg = _mapper.get_search_grid();  // search grid
  if (sg.stepping()[direction] > 0)
    return sg.origin()[direction];
  else return sg.origin()[direction] + sg.stepping()[direction] * sg.dims()[direction];
}

double GridIntersectionSearcher::get_highest_bound_(const size_t direction) const
{
  const auto & sg = _mapper.get_search_grid();  // search grid
  if (sg.stepping()[direction] < 0)
    return sg.origin()[direction];
  else return sg.origin()[direction] + sg.stepping()[direction] * sg.dims()[direction];
}


std::vector<size_t>
GridIntersectionSearcher::collision(const angem::LineSegment<double> & segment)
{
  const auto & sg = _mapper.get_search_grid();
  std::unordered_set<size_t> result;
  angem::CollisionGJK<double> checker;
  if (checker.check(segment, sg.get_bounding_box()))  // segment intersects the grid
  {
    // find first search location
    const auto start_end = find_first_and_last_location_(segment);
    const angem::Line<3,double> line(segment.line());

    auto cur = start_end.first;
    size_t search_idx = sg.find_cell(cur);
    std::unordered_set<size_t> processed;
    // line x = x0 + a * t;
    double t = get_t_(segment, cur);
    const double t_stop = std::min(get_t_(segment, segment.second()),
                                   get_t_(segment, start_end.second));
    int nit = 0;
    while (t <= t_stop)
    {
      const size_t search_cell = sg.find_cell(cur);
      process_(search_cell,result);
      t = step_(search_cell, cur, segment);
      cur =  segment.first() + line.direction() * t;
      nit++;
    }
  }
  return std::vector<size_t>(result.begin(), result.end());
}

std::pair<angem::Point<3,double>,angem::Point<3,double>>
GridIntersectionSearcher::
find_first_and_last_location_(const angem::LineSegment<double> & segment) const
{
  const auto & sg = _mapper.get_search_grid();
  std::vector<angem::Point<3,double>> intersection;
  angem::collision(segment, sg.get_bounding_box(), intersection);
  std::sort(intersection.begin(), intersection.end(),
  // auto loc = *std::min_element(intersection.begin(), intersection.end(),
                          [segment](const angem::Point<3,double> & p1, const angem::Point<3,double> & p2)
                          {
                            // return segment.first().distance(p1) < segment.second().distance(p2);
                            return segment.first().distance(p1) < segment.first().distance(p2);
                          });
  if (!sg.in_bounds(intersection[0]))  // might be slightly off due to roundoff
    intersection[0] += 1e-3 * segment.line().direction() * sg.stepping(0);
  if (!sg.in_bounds(intersection[1]))  // might be slightly off due to roundoff
    intersection[1] -= 1e-3 * segment.line().direction() * sg.stepping(0);
  assert(sg.in_bounds(intersection[0]));

  return {intersection[0], intersection[1]};
}

double GridIntersectionSearcher::get_t_(const angem::LineSegment<double> & segment, const angem::Point<3,double>& p) const
{
  const auto tangent = segment.line().direction();
  for (size_t i = 0; i < 3; ++i)
    if (std::fabs(tangent(i)) > 1e-10)
      return (p(i) - segment.first()(i)) / tangent(i);
  throw std::invalid_argument("should not be here");
}

void GridIntersectionSearcher::process_(const size_t search_cell,
                                        std::unordered_set<size_t> & result)
{
  for (const size_t cell_index : _mapper.mapping(search_cell) )
    if (_grid.cell(cell_index).is_active())
    {
      result.insert(cell_index);
    }
}

double GridIntersectionSearcher::step_(size_t search_cell,
                                       const angem::Point<3,double> & pos,
                                       const angem::LineSegment<double> & segment) const
{
  const double t = get_t_(segment, pos);
  const auto & sg = _mapper.get_search_grid();
  const auto ijk = sg.get_ijk(search_cell);
  const double t_stop = get_t_(segment, segment.second());
  const double eps = 1e-4;
  double t_new = t_stop + sg.stepping(0)*eps;
  const auto tangent = segment.line().direction();

  for (size_t i = 0; i < 3; ++i)
    if (std::isnormal(tangent(i)))  // nonzero
  {
    std::array<int, 3> candidate = ijk;
    if (sg.stepping(i) * tangent(i) > 0)  // same direction
      candidate[i] += 1;
    else candidate[i] -= 1;
    if (sg.is_valid_vertex(candidate[0], candidate[1], candidate[2]))
    {
      auto v = sg.vertex(candidate[0], candidate[1], candidate[2]);
      const double t_trial = (v[i] + sg.stepping(i)*eps - segment.first()[i]) / tangent(i);
      t_new = std::min(t_new, t_trial);
    }
  }
  return t_new;
}

std::vector<size_t> GridIntersectionSearcher::collision(const angem::Polygon<double> & polygon)
{
  const angem::Point<3,double> vertical(0,0,1);
  if (std::fabs(polygon.normal().dot(vertical) - 1.0) < 1e-6)
    throw std::invalid_argument("Horizontal fractures not supported");

  const auto & sg = _mapper.get_search_grid();
  double bot = polygon.support({0, 0, -1})[2];
  if (bot > this->top()) throw std::invalid_argument("fracture not in bounds");
  bot = std::max(bot, bottom() + 1e-2 * std::fabs(sg.stepping(2)));

  double top = polygon.support(vertical)[2];
  if (top < bottom()) throw std::invalid_argument("fracture not in bounds");
  top = std::min(top, this->top() - 1e-2 * std::fabs(sg.stepping(2)));

  const int start_k = std::max<int>((bot - bottom()) / sg.stepping(2), int(0));
  const int end_k = std::min<int>((top - this->top()) / sg.stepping(2), sg.nz()-1);
  std::unordered_set<size_t> result;
  for (size_t k = start_k; k <= end_k; ++k)
  {
    // z-center of the layer
    const double z = sg.origin()[2] + sg.stepping(2) * (k + 0.5);
    const angem::Plane<double> section_plane(/*origin=*/{0, 0, z}, /*normal*/vertical);
    std::vector<angem::Point<3,double>> intersection;
    angem::collision(polygon, section_plane, intersection);
    if (intersection.size() != 2)
    {
      std::cout << "z = " << z << std::endl;
      std::cout << "poly" << std::endl;
      for (auto p : polygon.get_points())
        std::cout << p << std::endl;
      std::cout << "section:" << std::endl;
      for (auto p : intersection)
        std::cout << p << std::endl;
    }
    assert(intersection.size() == 2);
    angem::LineSegment<double> segment(intersection[0], intersection[1]);
    for (const size_t icell : collision(segment))
      result.insert(icell);
  }
  return std::vector<size_t>(result.begin(), result.end());
}

}  // end namespace gprs_data
