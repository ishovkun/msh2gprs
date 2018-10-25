#pragma once

#include "SurfaceMesh.hpp"

namespace mesh
{

template<typename Scalar>
SurfaceMesh<Scalar> make_surface_mesh(const angem::Point<3,Scalar> & v1,
                                      const angem::Point<3,Scalar> & v2,
                                      const angem::Point<3,Scalar> & origin,
                                      const std::size_t n1,
                                      const std::size_t n2)
{
  std::vector<angem::Point<3,Scalar>> poly_points(4);
  angem::Point<3,Scalar> current_layer = {0, 0, 0};

  SurfaceMesh<Scalar> mesh;
  for (std::size_t j=0; j<n2; ++j)
  {
    poly_points[0] = origin + current_layer;
    poly_points[1] = origin + current_layer + v1 / n1;
    poly_points[2] = origin + current_layer + v1 / n1 + v2 / n2;
    poly_points[3] = origin + current_layer + v2 / n2;

    for (std::size_t i=0; i<n1; ++i)
    {
      angem::Polygon<Scalar> poly(poly_points);
      mesh.insert(poly);

      for (auto & p : poly_points)
        p += v1 / n1;
    }

    current_layer += v2 / n2;
  }

  return std::move(mesh);
}

}
