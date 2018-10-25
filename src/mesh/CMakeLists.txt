# mesh/CMakeLists.txt

ADD_LIBRARY(mesh
  SurfaceMesh.hpp
  utils.hpp
)

SET_TARGET_PROPERTIES (
    mesh
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_LINK_LIBRARIES(mesh angem)
