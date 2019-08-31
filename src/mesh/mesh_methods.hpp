#pragma once

#include "angem/Point.hpp"
#include "angem/PointSet.hpp"
#include "angem/Polyhedron.hpp"
#include "SurfaceMesh.hpp"
#include "Face.hpp"

#include <unordered_set>

#ifdef USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include <uint256/uint256_t.h>
#endif

namespace mesh
{

using Point = angem::Point<3,double>;

#ifdef USE_BOOST
using hash_type = boost::multiprecision::uint256_t;
#else
using hash_type = uint256_t;
#endif

using FaceMap = std::unordered_map<hash_type, Face>;


extern const std::size_t MAX_HASHED_VERTICES;
extern const int INTERNAL_FACE_ID;

std::size_t estimate_max_vertices(const int n_polygon_vertices);

Point get_element_center(const angem::PointSet<3,double> & vertices,
                         const std::vector<std::size_t>  & ivertices);

std::vector<Point> get_vertex_coordinates(const std::vector<angem::Point<3,double>> & vertices,
                                          const std::vector<std::size_t>            & ivertices);

std::vector<Point> get_vertex_coordinates(const angem::PointSet<3,double> * p_vertices,
                                          const std::vector<std::size_t>  & ivertices);

hash_type hash_value(const std::vector<std::size_t> & ivertices);

std::vector<std::size_t> invert_hash(const hash_type & hash);

int get_face_marker(const hash_type & hash,
                    const std::unordered_map<hash_type, int> & map_physical_faces);

std::vector<std::size_t> &
get_face_neighbors(const std::vector<std::size_t> & face,
                   std::unordered_map<hash_type, std::vector<std::size_t>> map_faces);

std::vector<std::size_t> face_vertices(const hash_type hash,
                                       angem::PointSet<3,double> & vertices);


std::vector<std::vector<std::size_t>>
get_face_indices(const angem::Polyhedron<double> & poly,
                 const angem::PointSet<3,double> & vertices);


std::unordered_set<std::size_t> find_internal_vertices(const SurfaceMesh<double> & msh);

}
