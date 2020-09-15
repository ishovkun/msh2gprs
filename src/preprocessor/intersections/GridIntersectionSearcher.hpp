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

 private:
  double get_lowest_bound_(const size_t direction) const;
  double get_highest_bound_(const size_t direction) const;
  double get_t_(const angem::LineSegment<double> & segment, const angem::Point<3,double>& p) const;
  angem::Point<3,double> find_first_location_(const angem::LineSegment<double> & segment) const;
  void process_(const size_t search_cell, std::unordered_set<size_t> & result);
  double step_(size_t search_cell,
               const angem::Point<3,double> & pos,
               const angem::LineSegment<double> & segment) const;

  const mesh::Mesh & _grid;
  Mapper _mapper;
};

}  // end namespace gprs_data
