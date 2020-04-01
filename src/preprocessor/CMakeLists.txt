# preprocessor/CMakeLists.txt

ADD_LIBRARY(preprocessor
  Preprocessor.cpp
  parsers/YamlParser.cpp
  gmsh_interface/GmshInterface.cpp
  gmsh_interface/FeValues.cpp
  gmsh_interface/element_types.cpp
  utils.cpp
  CellPropertyManager.cpp
  EmbeddedFractureManager.cpp
  DiscreteFractureManager.cpp
  BoundaryConditionManager.cpp
  WellManager.cpp
  Well.cpp
  discretization/flow/DiscretizationBase.cpp
  discretization/flow/DiscretizationTPFA.cpp
  discretization/flow/DiscretizationDFM.cpp
  discretization/flow/DiscretizationEDFM.cpp
  DoFManager.cpp
  OutputDataGPRS.cpp
  OutputDataVTK.cpp
  OutputDataPostprocessor.cpp
  VTKWriter.cpp
  MultiScaleDataMSRSB.cpp
  MultiScaleDataMech.cpp
  GridEntityNumberingManager.cpp
  UnionFind.cpp
  discretization/mechanics/DiscretizationDFEM.cpp
  discretization/mechanics/DFEMElement.cpp
  discretization/mechanics/JacobiPreconditioner.cpp
  discretization/mechanics/PolyhedralElementDirect.cpp
  discretization/mechanics/FEMFaceDoFManager.cpp
)

SET_TARGET_PROPERTIES (
    preprocessor
    PROPERTIES LINKER_LANGUAGE CXX
)

set(preprocessor_libs ${preprocessor_libs} angem muparser mesh yaml-cpp)

if(USE_GMSH)
  TARGET_INCLUDE_DIRECTORIES(preprocessor PRIVATE ${GMSH_INCLUDE_DIRS})
  set(preprocessor_libs ${preprocessor_libs} ${Gmsh_LIBRARIES})
endif()

if(USE_EIGEN)
  TARGET_INCLUDE_DIRECTORIES(preprocessor PRIVATE ${EIGEN_INCLUDE_DIR})
  set(preprocessor_libs ${preprocessor_libs} Eigen3::Eigen)
endif()

if(METIS_FOUND)
  set(preprocessor_libs ${preprocessor_libs} ${METIS_LIBRARIES})
endif()

message(" linking preproc libs ${preprocessor_libs}")
TARGET_LINK_LIBRARIES(preprocessor ${preprocessor_libs})
