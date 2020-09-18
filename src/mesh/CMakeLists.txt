# mesh/CMakeLists.txt

ADD_LIBRARY(mesh
  ShapeID.cpp
  Mesh.cpp
  Cell.cpp
  Face.cpp
  # face_iterator.cpp
  # face_const_iterator.cpp
  # mesh_methods.cpp
  surface_mesh_methods.cpp
  active_cell_iterator.cpp
  active_cell_const_iterator.cpp
  active_face_const_iterator.cpp
  Subdivision.cpp
  FaceSplitter.cpp
  CellSplitter.cpp
  FractureTip.cpp
  MeshStatsComputer.cpp
  CartesianMeshBuilder.cpp
  RefinementAspectRatio.cpp
  io/VTKReader.cpp
  io/VTKWriter.cpp
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
