#pragma once

#include "Mesh.hpp"

namespace mesh {

/*
** This class figures out whether fracture vertices belong to the tip or not
 */
class FractureTip
{
 public:
  FractureTip(const SurfaceMesh<double> & face_mesh,
              const Mesh & grid,
              std::unordered_map<size_t,size_t> & _surface_vertex_to_global);

  virtual ~FractureTip() = default;

  bool contains(const size_t vertex);

 private:
  std::unordered_set<size_t> find_global_boundary_();

  const SurfaceMesh<double> & _face_mesh;
  const Mesh & _grid;
  std::unordered_map<size_t,size_t> & _surface_vertex_to_global;
  std::vector<std::vector<size_t>> _tip_verts;
  std::unordered_set<size_t> _global_boundary;
};

}  // end namespace mesh
