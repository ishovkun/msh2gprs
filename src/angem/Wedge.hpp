#pragma once

#include <Polyhedron.hpp>

namespace angem
{

/*
  Node numbering:
  Prism:                      Prism15:               Prism18:

            w
            ^
            |
            3                       3                      3
          ,/|`\                   ,/|`\                  ,/|`\
        ,/  |  `\               12  |  13              12  |  13
      ,/    |    `\           ,/    |    `\          ,/    |    `\
      4------+------5         4------14-----5        4------14-----5
      |      |      |         |      8      |        |      8      |
      |    ,/|`\    |         |      |      |        |    ,/|`\    |
      |  ,/  |  `\  |         |      |      |        |  15  |  16  |
      |,/    |    `\|         |      |      |        |,/    |    `\|
    ,|      |      |\        10     |      11       10-----17-----11
  ,/ |      0      | `\      |      0      |        |      0      |
  u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
      |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
      |,/         `\|         |,/         `\|        |,/         `\|
      1-------------2         1------9------2        1------9------2

*/

template<typename Scalar>
class Wedge: public Polyhedron<Scalar>
{
 public:
  // CONSTRICTORS
  Wedge();
  // input given with vtk numbering
  Wedge(const std::vector<Point<3,Scalar>> & vertices,
        const std::vector<std::size_t>     & indices);
  // SETTERS
  virtual void set_data(const std::vector<Point<3,Scalar>> & vertices) override;
  void set_data(const std::vector<Point<3,Scalar>> & vertices,
                const std::vector<std::size_t>     & indices);
};


template<typename Scalar>
Wedge<Scalar>::Wedge()
    :
    Polyhedron<Scalar>(13)
{}


template<typename Scalar>
Wedge<Scalar>::Wedge(const std::vector<Point<3,Scalar>> & vertices,
                     const std::vector<std::size_t>     & indices)
    :
    Polyhedron<Scalar>(13)
{
  set_data(vertices, indices);
}


template<typename Scalar>
void
Wedge<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices)
{
  if (vertices.size() == 15 or vertices.size() == 18)
    throw NotImplemented("Only first order Wedges implemented");
  assert(vertices.size() == 6);

  this->points = vertices;
  this->faces.resize(5);
  this->faces[0] = {0, 1, 2};
  this->faces[1] = {3, 4, 5};
  this->faces[2] = {0, 3, 4, 1};
  this->faces[3] = {1, 2, 5, 4};
  this->faces[4] = {0, 3, 5, 2};
}


template<typename Scalar>
void
Wedge<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices,
                        const std::vector<std::size_t>     & indices)
{
  std::vector<Point<3,Scalar>> points(indices.size());
  int counter = 0;
  for (const auto & ivert : indices)
  {
    points[counter] = vertices[ivert];
    counter++;
  }

  set_data(points);
}

}
