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
  // GETTERS
  virtual Scalar volume() const override;
  static Scalar volume(const Point<3,Scalar> & v1,
                       const Point<3,Scalar> & v2,
                       const Point<3,Scalar> & v3,
                       const Point<3,Scalar> & v4);

  // don't create polyhedron, just give me vector of faces with indices in
  // the global vector
  static std::vector<std::vector<std::size_t>>
  get_faces(const std::vector<std::size_t>  & indices);
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


template<typename Scalar>
std::vector<std::vector<std::size_t>>
Tetrahedron<Scalar>::get_faces(const std::vector<std::size_t> & indices)
{
  std::vector<std::vector<std::size_t>> global_faces(4);
  global_faces[0] = {indices[0], indices[1], indices[3]};
  global_faces[1] = {indices[1], indices[2], indices[3]};
  global_faces[2] = {indices[0], indices[2], indices[3]};
  global_faces[3] = {indices[0], indices[1], indices[2]};
  return std::move(global_faces);
}


template<typename Scalar>
Scalar Tetrahedron<Scalar>::volume() const
{
  // volume = 1 / 6 |a·(b x c)|
  Point<3,Scalar> a = this->points[1] - this->points[0];
  Point<3,Scalar> b = this->points[2] - this->points[0];
  Point<3,Scalar> c = this->points[3] - this->points[0];
  return static_cast<Scalar>(1./6.) * a.dot(b.cross(c));
}


template<typename Scalar>
inline
Scalar Tetrahedron<Scalar>::volume(const Point<3,Scalar> & v0,
                                   const Point<3,Scalar> & v1,
                                   const Point<3,Scalar> & v2,
                                   const Point<3,Scalar> & v3)
{
  // volume = 1 / 6 |a·(b x c)|
  Point<3,Scalar> a = v1 - v0;
  Point<3,Scalar> b = v2 - v0;
  Point<3,Scalar> c = v3 - v0;
  return static_cast<Scalar>(1./6.) * fabs(a.dot(b.cross(c)));
}

}
