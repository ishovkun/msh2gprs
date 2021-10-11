#pragma once
#include "mesh/Mesh.hpp"
#include "Mapper.hpp"
#include "angem/LineSegment.hpp"

namespace gprs_data {

class GridIntersectionSearcher {
 public:
  GridIntersectionSearcher(const mesh::Mesh & grid);

  double top() const noexcept;
  double bottom() const noexcept;
  // find collision of a line segment with the grid (for wells)
  // Input:
  // [p1, p2] : defines the line segment
  // returns:
  // list of cell indices intersected by the segment
  std::vector<size_t> collision(const angem::LineSegment<double> & segment);
  // find collision of a polygon with the grid (for fractures)
  // Input:
  // [p1, p2] : defines the line segment
  // returns:
  // list of cell indices intersected by the segment
  std::vector<size_t> collision(const angem::Polygon<double> & polygon);
  // find cell that contains the point
  // Input:
  // \param[in] p : defines the point
  // Returns the index of the cell that contains the point. If the point is
  // out of the domain, returns _grid.n_vertices()
  size_t find_cell(angem::Point<3,double> const & p) const;

 private:
  double get_lowest_bound_(const size_t direction) const;
  double get_highest_bound_(const size_t direction) const;
  double get_t_(const angem::LineSegment<double> & segment, const angem::Point<3,double>& p) const;
  std::pair<angem::Point<3,double>,angem::Point<3,double>>
  find_first_and_last_location_(const angem::LineSegment<double> & segment) const;
  void process_(const size_t search_cell, std::unordered_set<size_t> & result);
  double step_(size_t search_cell,
               const angem::Point<3,double> & pos,
               const angem::LineSegment<double> & segment) const;

  const mesh::Mesh & _grid;
  Mapper _mapper;
};

}  // end namespace gprs_data
