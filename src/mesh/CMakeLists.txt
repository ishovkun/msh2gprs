# mesh/CMakeLists.txt

ADD_LIBRARY(mesh
  SurfaceMesh.hpp
  utils.hpp
  Mesh.hpp
  ShapeID.hpp
  cell_iterator.hpp
  face_iterator.hpp
  mesh_methods.hpp
  surface_mesh_methods.hpp
  # implementation
  ShapeID.cpp
  Mesh.cpp
  cell_iterator.cpp
  face_iterator.cpp
  mesh_methods.cpp
  surface_mesh_methods.cpp
)

SET_TARGET_PROPERTIES (
    mesh
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC
	${CMAKE_SOURCE_DIR}/src
)


# # TARGET_LINK_LIBRARIES(mesh angem uint256)
# TARGET_LINK_LIBRARIES(mesh angem vli)
TARGET_LINK_LIBRARIES(mesh angem ${Boost_LIBRARIES})
