# mesh/CMakeLists.txt

ADD_LIBRARY(mesh
  ShapeID.cpp
  Mesh.cpp
  cell_iterator.cpp
  const_cell_iterator.cpp
  face_iterator.cpp
  const_face_iterator.cpp
  mesh_methods.cpp
  surface_mesh_methods.cpp
)

SET_TARGET_PROPERTIES (
    mesh
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/uint256
)


if(Boost_FOUND)
  TARGET_LINK_LIBRARIES(mesh angem ${Boost_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(mesh angem uint256)
endif()
