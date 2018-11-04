# mesh/CMakeLists.txt

ADD_LIBRARY(mesh
  SurfaceMesh.hpp
  utils.hpp
  Mesh.hpp
  # implementation
  Mesh.cpp
)

SET_TARGET_PROPERTIES (
    mesh
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC
	${CMAKE_SOURCE_DIR}/src
)


TARGET_LINK_LIBRARIES(mesh angem uint256)
