#pragma once

#include <Polyhedron.hpp>
#include <Tetrahedron.hpp>
#include <Exceptions.hpp>

/* Hexahedron vertex numbering:
Hexahedron:             Hexahedron20:          Hexahedron27:

       v
3----------2            3----13----2           3----13----2
|\     ^   |\           |\         |\          |\         |\
| \    |   | \          | 15       | 14        |15    24  | 14
|  \   |   |  \         9  \       11 \        9  \ 20    11 \
|   7------+---6        |   7----19+---6       |   7----19+---6
|   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
 \  |    \  \  |         \  17      \  18       \ 17    25 \  18
  \ |     \  \ |         10 |        12|        10 |  21    12|
   \|      w  \|           \|         \|          \|         \|
    4----------5            4----16----5           4----16----5
*/

namespace angem
{

template<typename Scalar>
class Hexahedron: public Polyhedron<Scalar>
{
  public:
  // CONSTRICTORS
  Hexahedron();
  Hexahedron(const std::vector<Point<3,Scalar>> & vertices,
             const std::vector<std::size_t>     & indices);
  // SETTERS
  virtual void set_data(const std::vector<Point<3,Scalar>> & vertices) override;
  void set_data(const std::vector<Point<3,Scalar>> & vertices,
                const std::vector<std::size_t>     & indices);
  // GETTERS
  virtual Scalar volume() const override;
  // don't create polyhedron, just give me vector of faces with indices in
  // the global vector
  static std::vector<std::vector<std::size_t>>
  get_faces(const std::vector<std::size_t>     & indices);
};


template<typename Scalar>
Hexahedron<Scalar>::Hexahedron()
    :
    Polyhedron<Scalar>(12)
{}


template<typename Scalar>
Hexahedron<Scalar>::Hexahedron(const std::vector<Point<3,Scalar>> & vertices,
                               const std::vector<std::size_t>     & indices)
    :
    Polyhedron<Scalar>(12)
{
  set_data(vertices, indices);
}


template<typename Scalar>
void
Hexahedron<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices)
{
  if (vertices.size() == 20 or vertices.size() == 27)
    throw NotImplemented("Only first order Hexas implemented");
  assert(vertices.size() == 8);

  this->points = vertices;
  this->faces.resize(6);
  this->faces[0] = {0, 1, 2, 3};
  this->faces[1] = {4, 5, 6, 7};
  this->faces[2] = {0, 4, 5, 1};
  this->faces[3] = {3, 7, 6, 2};
  this->faces[4] = {0, 4, 7, 3};
  this->faces[5] = {1, 5, 6, 2};
}


template<typename Scalar>
void
Hexahedron<Scalar>::set_data(const std::vector<Point<3,Scalar>> & vertices,
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


template<typename Scalar>
std::vector<std::vector<std::size_t>>
Hexahedron<Scalar>::get_faces(const std::vector<std::size_t>     & indices)
{
  std::vector<std::vector<std::size_t>> global_faces(6);
  global_faces[0] = {indices[0], indices[1], indices[2], indices[3]};
  global_faces[1] = {indices[4], indices[5], indices[6], indices[7]};
  global_faces[2] = {indices[0], indices[4], indices[5], indices[1]};
  global_faces[3] = {indices[3], indices[7], indices[6], indices[2]};
  global_faces[4] = {indices[0], indices[4], indices[7], indices[3]};
  global_faces[5] = {indices[1], indices[5], indices[6], indices[2]};
  return std::move(global_faces);
}


template<typename Scalar>
Scalar Hexahedron<Scalar>::volume() const
{
  /* volume = sum of volumes of 6 tetrahedrons
   * see
   * See discussions, stats, and author profiles for this publication at:
   * https://www.researchgate.net/publication/260322098
   *Reentry Flows in Chemical Non-Equilibrium in Three-Dimensions */

  const auto & points = this->points;
  return Tetrahedron<Scalar>::volume(points[3], points[6], points[0], points[4]) +
      Tetrahedron<Scalar>::volume(points[4], points[0], points[5], points[6]) +
      Tetrahedron<Scalar>::volume(points[0], points[1], points[6], points[5]) +
      Tetrahedron<Scalar>::volume(points[0], points[1], points[6], points[2]) +
      Tetrahedron<Scalar>::volume(points[3], points[2], points[6], points[0]) +
      Tetrahedron<Scalar>::volume(points[7], points[3], points[6], points[0]);
}

}
