#pragma once

#include "Point.hpp"
#include "Plane.hpp"

namespace angem
{

// generic template
template <int dim, typename Scalar>
class Polygon
{
 public:
  Polygon();
  Polygon(std::vector< Point<dim,Scalar> > & bounding_points);
  void set_data(std::vector< Point<dim,Scalar> > & bounding_points);

 protected:
  std::vector< Point<dim,Scalar> > bounding_points;
};

// ----------------------------- 2D ------------------------------------------
template <typename Scalar>
class Polygon<2,Scalar>
{
 public:
  Polygon();
  Polygon(std::vector< Point<2,Scalar> > & bounding_points);
  void set_data(std::vector< Point<2,Scalar> > & bounding_points);

 protected:
  std::vector< Point<2,Scalar> > bounding_points;
};


template <typename Scalar>
Polygon<2,Scalar>::Polygon()
{}


template <typename Scalar>
Polygon<2,Scalar>::Polygon(std::vector< Point<2,Scalar> > & bounding_points)
{
  set_data(bounding_points);
}


template <typename Scalar>
void
Polygon<2,Scalar>::set_data(std::vector< Point<2,Scalar> > & v_bounding_points)
{
  assert(bounding_points.size() > 2);

  bounding_points = v_bounding_points;
}
// ----------------------------- 3D ------------------------------------------
template <typename Scalar>
class Polygon<3,Scalar>
{
 public:
  Polygon();
  Polygon(std::vector< Point<3,Scalar> > & bounding_points);
  void set_data(std::vector< Point<3,Scalar> > & bounding_points);
  Polygon<2,Scalar> to_2D() const;

  Plane<Scalar> plane;

 protected:
  std::vector< Point<3,Scalar> > bounding_points;
};


template <typename Scalar>
Polygon<3,Scalar>::Polygon()
{}


template <typename Scalar>
Polygon<3,Scalar>::Polygon(std::vector< Point<3,Scalar> > & bounding_points)
{
  set_data(bounding_points);
}


template <typename Scalar>
void Polygon<3,Scalar>::set_data(std::vector< Point<3,Scalar> > & v_bounding_points)
{
  assert(bounding_points.size() > 2);
  assert(align_on_plane(v_bounding_points));

  bounding_points = v_bounding_points;

  plane = Plane<Scalar>(bounding_points[0],
                        bounding_points[1],
                        bounding_points[2]);
}


template <typename Scalar>
Polygon<2,Scalar> Polygon<3,Scalar>::to_2D() const
{
  std::vector<Point<2,Scalar>> points(bounding_points.size());
  for (std::size_t i=0; i< bounding_points.size(); ++i)
  {
    const Point<3,Scalar> lc = plane.local_coordinates(bounding_points[i]);
    points[i] = Point<2,Scalar>(lc(0), lc(1));
  }
  return Polygon<2,Scalar>(points);
}

}  // end namespace
