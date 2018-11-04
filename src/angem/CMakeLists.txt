# angem/CMakeLists.txt



ADD_LIBRARY(angem
  Point.hpp
  PointSet.hpp
  PointSet.cpp
  Basis.hpp
  Plane.hpp
  Shape.hpp
  Polygon.hpp
  Polyhedron.hpp
  Rectangle.hpp
  Collisions.hpp
  CollisionGJK.hpp
  Line.hpp
  PolyGroup.hpp
)

SET_TARGET_PROPERTIES (
    angem
    PROPERTIES LINKER_LANGUAGE CXX
)
