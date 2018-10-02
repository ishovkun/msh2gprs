#include <Point.hpp>
#include <Polygon.hpp>

namespace angem
{

template <typename Scalar>
struct PolyGroup
{
  std::vector<Point<3,double>> vertices;
  std::vector<Polygon<Scalar>> polygons;
};

}
