#include <Point.hpp>
#include <Polygon.hpp>

namespace angem
{

template <typename Scalar>
struct PolyGroup
{
  std::vector<Point<3,double>> vertices;
  std::vector<Polygon<Scalar>> polygons;
  PolyGroup<Scalar> & operator=(const PolyGroup<Scalar> &other);
};


template <typename Scalar>
PolyGroup<Scalar> &
PolyGroup<Scalar>::operator=(const PolyGroup<Scalar> &other)
{
  vertices = other.vertices;

  for (const auto & poly : other.polygons)
  {
    std::vector<Point<3,Scalar>*> p_points;
    for (const auto p_p : poly.get_points())
    {
      const std::size_t ind = std::distance(other.vertices.begin(),
                                            std::find(other.vertices.begin(),
                                                      other.vertices.end(),
                                                      *p_p));
      p_points.push_back(&vertices[ind]);
    }
    polygons.push_back(Polygon<Scalar>(p_points));
  }
}


}
