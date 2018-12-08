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
  static std::unique_ptr<Polyhedron<Scalar>>
  create(const std::vector<Point<3,Scalar>> & vertices,
         const std::vector<std::size_t>     & indices,
         const int                            vtk_id = -1)
  {
    switch (vtk_id)
    {
      case 10:
        return std::make_unique<Tetrahedron<Scalar>>(vertices, indices);
        break;
      case 12:
        return std::make_unique<Hexahedron<Scalar>>(vertices, indices);
        break;
      case 13:
        return std::make_unique<Wedge<Scalar>>(vertices, indices);
        break;
      case 14:
        return std::make_unique<Pyramid<Scalar>>(vertices, indices);
        break;
      default:
        {
          // try to construct based on number of points
          const int n_verts = vertices.size();
          switch (n_verts)
          {
            case 4:
              return std::make_unique<Tetrahedron<Scalar>>(vertices, indices);
              break;
            case 5:
              return std::make_unique<Hexahedron<Scalar>>(vertices, indices);
              break;
            case 6:
              return std::make_unique<Wedge<Scalar>>(vertices, indices);
              break;
            case 8:
              return std::make_unique<Pyramid<Scalar>>(vertices, indices);
              break;
            default:
              throw NotImplemented("3D element does not exist");
          }
        }
    }
  }

  template<typename Scalar>
  static std::unique_ptr<Polyhedron<Scalar>>
  create(const int  vtk_id = -1)
  {
    switch (vtk_id)
    {
      case 10:
        return std::make_unique<Tetrahedron<Scalar>>();
      case 12:
        return std::make_unique<Hexahedron<Scalar>>();
      case 13:
        return std::make_unique<Wedge<Scalar>>();
      case 14:
        return std::make_unique<Pyramid<Scalar>>();
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
