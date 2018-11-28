#pragma once

// #include <Polyhedron.hpp>
#include <Tetrahedron.hpp>
#include <Hexahedron.hpp>
#include <Wedge.hpp>
#include <Pyramid.hpp>
#include <Exceptions.hpp>

namespace angem
{

struct NotImplemented;

class PolyhedronFactory
{
 public:
  template<typename Scalar>
  static Polyhedron<Scalar>
  create(const std::vector<Point<3,Scalar>> & vertices,
         const std::vector<std::size_t>     & indices,
         const int                            vtk_id = -1)
  {
    switch (vtk_id)
    {
      case 10:
        return Tetrahedron<Scalar>(vertices, indices);
      case 12:
        return Hexahedron<Scalar>(vertices, indices);
      case 13:
        return Wedge<Scalar>(vertices, indices);
      case 14:
        return Pyramid<Scalar>(vertices, indices);
      default:
        {
          // try to construct based on number of points
          const int n_verts = vertices.size();
          switch (n_verts)
          {
            case 4:
              return Tetrahedron<Scalar>(vertices, indices);
            case 5:
              return Pyramid<Scalar>(vertices, indices);
            case 6:
              return Wedge<Scalar>(vertices, indices);
            case 8:
              return Hexahedron<Scalar>(vertices, indices);
          }
          throw NotImplemented("3D element does not exist");
        }
    }
  }

  template<typename Scalar>
  static Polyhedron<Scalar>
  create(const int  vtk_id = -1)
  {
    switch (vtk_id)
    {
      case 10:
        return Tetrahedron<Scalar>();
      case 12:
        return Hexahedron<Scalar>();
      case 13:
        return Wedge<Scalar>();
      case 14:
        return Pyramid<Scalar>();
      default:
        throw NotImplemented("3D element does not exist");
    }
  }


  template<typename Scalar>
  static std::vector<std::vector<std::size_t>>
  get_global_faces(const std::vector<std::size_t> & indices,
                   const int  vtk_id)
  {
    switch (vtk_id)
    {
      case 10:
        return Tetrahedron<Scalar>::get_faces(indices);
      case 12:
        return Hexahedron<Scalar>::get_faces(indices);
      case 13:
        return Wedge<Scalar>::get_faces(indices);
      case 14:
        return Pyramid<Scalar>::get_faces(indices);
      default:
        throw NotImplemented("3D element does not exist");
    }
  }


 private:
  PolyhedronFactory(){};
};

}