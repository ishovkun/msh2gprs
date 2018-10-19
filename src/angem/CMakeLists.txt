# angem/CMakeLists.txt



ADD_LIBRARY(angem
  Point.hpp
  PointSet.hpp
  PointSet.cpp
  Basis.hpp
  Plane.hpp
  Shape.hpp
  Rectangle.hpp
  Collisions.hpp
  CollisionGJK.hpp
  Line.hpp
  PolyGroup.hpp
  SurfaceMesh.hpp
)

SET_TARGET_PROPERTIES (
    angem
    PROPERTIES LINKER_LANGUAGE CXX
)
