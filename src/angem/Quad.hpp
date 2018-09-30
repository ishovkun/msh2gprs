#pragma once
#include <Polygon.hpp>

namespace angem
{

template<typename Scalar>
class Quad: public Polygon<Scalar>
{
 public:
  Quad(const std::vector<Point<3,Scalar>> & points);
  const std::vector<Point<3,Scalar>> & get_points() {return v_points;}

 protected:
  std::vector<Point<3,Scalar>> v_points;
};


template<typename Scalar>
Quad<Scalar>:: Quad(const std::vector<Point<3,Scalar>> & points)
    :
    Polygon<Scalar>::Polygon()
{
  assert(points.size() == 4);
  // i want v_points to have a proper vtk numbering
  Point<3,Scalar> center = compute_center_mass(points);
  // std::cout << "center = " << center << std::endl;
  // for (auto & p : points)
  //   std::cout << p << std::endl;

  std::vector<Point<3,Scalar>> copy = points;
  v_points.push_back(copy.front());
  copy.erase(copy.begin());

  while (!copy.empty())
  {
    for (std::size_t i=0; i<copy.size(); ++i)
    {
      if (!(  (copy[i] - center).parallel(v_points.back() - center)  ))
      {
        v_points.push_back(copy[i]);
        copy.erase(copy.begin() + i);
        break;
      }
    }
  }

  Polygon<Scalar>::set_data(v_points);
}



}
