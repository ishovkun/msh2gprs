#include "Point.hpp"
#include "Polygon.hpp"

namespace angem
{

template <typename Scalar>
class Cell
{
 public:
  Cell(const std::vector<Point<3,Scalar>> & vertices);


  std::vector<Point<3,Scalar>> vertices;
}

}
