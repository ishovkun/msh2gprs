#include <Point.hpp>
#include <Polygon.hpp>

namespace angem
{

template <typename Scalar>
struct PolyGroup
{
  std::vector<Point<3,Scalar>>          vertices;
  std::vector<std::vector<std::size_t>> polygons;  // indices

  void add(const PolyGroup<Scalar> & other);
  void add(const Polygon<Scalar> & poly);
};


template <typename Scalar>
void PolyGroup<Scalar>::add(const PolyGroup<Scalar> & other)
{
  // map vertices in other to new vertices in this
  std::vector<std::size_t> indices;
  for (const auto & p : other.vertices)
  {
    const std::size_t ind = insert(p, vertices, 1e-6);
    indices.push_back(ind);
  }

  for (const auto & poly : other.polygons)
  {
    std::vector<std::size_t> poly_new;
    for (std::size_t i : poly)
      poly_new.push_back(indices[i]);

    polygons.push_back(std::move(poly_new));
  }
}


template <typename Scalar>
void PolyGroup<Scalar>::add(const Polygon<Scalar> & poly)
{
  for (const auto & p : poly.get_points())
  {
    std::vector<std::size_t> indices;
    const std::size_t ind = insert(p, vertices, 1e-6);
    indices.push_back(ind);
  }
}

}
