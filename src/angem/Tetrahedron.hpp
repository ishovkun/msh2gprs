#pragma once

#include <Polyhedron.hpp>
#include <Exceptions.hpp>

namespace angem
{

/*
Node numbering:
Tetrahedron:                          Tetrahedron10 (not implemented):

                   v
                 .
               ,/
              /
           2                                     2
         ,/|`\                                 ,/|`\
       ,/  |  `\                             ,/  |  `\
     ,/    '.   `\                         ,6    '.   `5
   ,/       |     `\                     ,/       8     `\
 ,/         |       `\                 ,/         |       `\
0-----------'.--------1 --> u         0--------4--'.--------1
 `\.         |      ,/                 `\.         |      ,/
    `\.      |    ,/                      `\.      |    ,9
       `\.   '. ,/                           `7.   '. ,/
          `\. |/                                `\. |/
             `3                                    `3
                `\.
                   ` w
 */

template<typename Scalar>
class Tetrahedron: public Polyhedron<Scalar>
{
 public:
  // CONSTRICTORS
  Tetrahedron();
  // input given with vtk numbering
  Tetrahedron(const std::vector<Point<3,Scalar>> & vertices,
              const std::vector<std::size_t>     & indices);
  // SETTERS
  virtual void set_data(const std::vector<Point<3,Scalar>> & vertices) override;
  void set_data(const std::vector<Point<3,Scalar>> & vertices,
                const std::vector<std::size_t>     & indices);
};


template<typename Scalar>
Tetrahedron<Scalar>::Tetrahedron()
    :
    Polyhedron<Scalar>(10)
{}


template<typename Scalar>
Tetrahedron<Scalar>::Tetrahedron(const std::vector<Point<3,Scalar>> & vertices,
                                 const std::vector<std::size_t>     & indices)
    :
    Polyhedron<Scalar>(10)
{
  set_data(vertices, indices);
}


template<typename Scalar>
void
Tetrahedron<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices)
{
  if (vertices.size() == 10)
    throw NotImplemented("Only first order Tetras implemented");
  assert(vertices.size() == 4);

  this->points = vertices;
  this->faces.resize(4);
  this->faces[0] = {0, 1, 3};
  this->faces[1] = {1, 2, 3};
  this->faces[2] = {0, 2, 3};
  this->faces[3] = {0, 1, 2};
}


template<typename Scalar>
void
Tetrahedron<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices,
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
