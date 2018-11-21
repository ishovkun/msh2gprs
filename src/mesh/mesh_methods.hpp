#pragma once

#include <angem/Point.hpp>
#include <angem/PointSet.hpp>
#include <angem/Polyhedron.hpp>
#include <SurfaceMesh.hpp>
#include <uint256/uint256_t.h>

#include <unordered_set>


namespace mesh
{

using Point = angem::Point<3,double>;

extern const std::size_t MAX_HASHED_VERTICES;
extern const int INTERNAL_FACE_ID;

std::size_t estimate_max_vertices(const int n_polygon_vertices);

Point get_element_center(const angem::PointSet<3,double> & vertices,
                         const std::vector<std::size_t>  & ivertices);

std::vector<Point> get_vertex_coordinates(const angem::PointSet<3,double> & vertices,
                                          const std::vector<std::size_t>  & ivertices);

uint256_t hash_value(const std::vector<std::size_t> & ivertices);

std::vector<std::size_t> invert_hash(const uint256_t & hash);

int get_face_marker(const uint256_t & hash,
                    const std::unordered_map<uint256_t, int> & map_physical_faces);

std::vector<std::size_t> &
get_face_neighbors(const std::vector<std::size_t> & face,
                   std::unordered_map<uint256_t, std::vector<std::size_t>> map_faces);

std::vector<std::size_t> face_vertices(const uint256_t hash,
                                       angem::PointSet<3,double> & vertices);


std::vector<std::vector<std::size_t>>
get_face_indices(const angem::Polyhedron<double> & poly,
                 const angem::PointSet<3,double> & vertices);


std::unordered_set<std::size_t> find_internal_vertices(const SurfaceMesh<double> & msh);

}
