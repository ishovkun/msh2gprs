#include <Point.hpp>
#include <Polygon.hpp>
#include <PointSet.hpp>

namespace angem
{

template <typename Scalar>
struct PolyGroup
{
  // std::vector<Point<3,Scalar>>          vertices;
  PointSet<3,Scalar>                    vertices;
  std::vector<std::vector<std::size_t>> polygons;  // indices
  std::vector<int>                      markers;

  void add(const PolyGroup<Scalar> & other);
  void add(const Polygon<Scalar> & poly,
           const int marker = 0);
};


template <typename Scalar>
void PolyGroup<Scalar>::add(const PolyGroup<Scalar> & other)
{
  assert(other.markers.size() == other.polygons.size());

  // map vertices in other to new vertices in this
  std::vector<std::size_t> indices;
  for (const auto & p : other.vertices)
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }

  for (const auto & poly : other.polygons)
  {
    std::vector<std::size_t> poly_new;
    for (std::size_t i : poly)
      poly_new.push_back(indices[i]);

    polygons.push_back(std::move(poly_new));
  }

  for (const int marker : other.markers)
    markers.push_back(marker);
}


template <typename Scalar>
void PolyGroup<Scalar>::add(const Polygon<Scalar> & poly,
                            const int               marker)
{
  std::vector<std::size_t> indices;
  for (const auto & p : poly.get_points())
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }
  polygons.push_back(indices);

  markers.push_back(marker);
}

}
