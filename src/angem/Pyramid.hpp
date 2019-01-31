#pragma once

#include <Polyhedron.hpp>

/*
Pyramid:                     Pyramid13:                   Pyramid14:

               4                            4                            4
             ,/|\                         ,/|\                         ,/|\
           ,/ .'|\                      ,/ .'|\                      ,/ .'|\
         ,/   | | \                   ,/   | | \                   ,/   | | \
       ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
     ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
   ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
 ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
 `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
   `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
     `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \
       `\.'      `\     `\`          `\.'             `\`         `\.'             `\`
          1----------------2            1--------8-------2           1--------8-------2
                    `\
                       u
*/

namespace angem
{

template<typename Scalar>
class Pyramid : public Polyhedron<Scalar>
{
 public:
  // CONSTRICTORS
  Pyramid();
  // input given with vtk numbering
  Pyramid(const std::vector<Point<3,Scalar>> & vertices,
          const std::vector<std::size_t>     & indices);
  // SETTERS
  void set_data(const std::vector<Point<3,Scalar>> & vertices);
  void set_data(const std::vector<Point<3,Scalar>> & vertices,
                const std::vector<std::size_t>     & indices);

  // don't create polyhedron, just give me vector of faces with indices in
  // the global vector
  static std::vector<std::vector<std::size_t>>
  get_faces(const std::vector<std::size_t>     & indices);
};


template<typename Scalar>
Pyramid<Scalar>::Pyramid()
    :
    Polyhedron<Scalar>(14)
{}


template<typename Scalar>
Pyramid<Scalar>::Pyramid(const std::vector<Point<3,Scalar>> & vertices,
                         const std::vector<std::size_t>     & indices)
    :
    Polyhedron<Scalar>(14)
{
  set_data(vertices, indices);
}


template<typename Scalar>
void
Pyramid<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices)
{
  if (vertices.size() == 13 or vertices.size() == 14)
    throw NotImplemented("Only first order Pyramids implemented");
  assert(vertices.size() == 5);

  this->points = vertices;
  this->faces.resize(5);
  this->faces[0] = {0, 1, 2, 3};
  this->faces[1] = {0, 4, 3};
  this->faces[2] = {0, 1, 4};
  this->faces[3] = {1, 4, 2};
  this->faces[4] = {2, 3, 4};
}


template<typename Scalar>
void
Pyramid<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices,
                          const std::vector<std::size_t>     & indices)
{
  if (vertices.size() == 13 or vertices.size() == 14)
    throw NotImplemented("Only first order Pyramids implemented");
  assert(indices.size() == 5);

  this->points.resize(indices.size());
  int counter = 0;
  for (const auto & ivert : indices)
  {
    this->points[counter] = vertices[ivert];
    counter++;
  }

  this->faces.resize(5);
  this->faces[0] = {0, 1, 2, 3};
  this->faces[1] = {0, 4, 3};
  this->faces[2] = {0, 1, 4};
  this->faces[3] = {1, 4, 2};
  this->faces[4] = {2, 3, 4};
}


template<typename Scalar>
std::vector<std::vector<std::size_t>>
Pyramid<Scalar>::get_faces(const std::vector<std::size_t>     & indices)
{
  std::vector<std::vector<std::size_t>> global_faces(5);
  global_faces[0] = {indices[0], indices[1], indices[2], indices[3]};
  global_faces[1] = {indices[0], indices[4], indices[3]};
  global_faces[2] = {indices[0], indices[1], indices[4]};
  global_faces[3] = {indices[1], indices[4], indices[2]};
  global_faces[4] = {indices[2], indices[3], indices[4]};
  return global_faces;
}

}
