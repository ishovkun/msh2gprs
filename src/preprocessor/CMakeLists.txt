# preprocessor/CMakeLists.txt

ADD_LIBRARY(preprocessor
  Preprocessor.cpp
  parsers/YamlParser.cpp
  parsers/GmshReader.cpp
  CellPropertyManager.cpp
  EmbeddedFractureManager.cpp
  DiscreteFractureManager.cpp
  BoundaryConditionManager.cpp
  WellManager.cpp
  Well.cpp
  discretization/DiscretizationBase.cpp
  discretization/DiscretizationTPFA.cpp
  discretization/DiscretizationDFM.cpp
  discretization/DiscretizationEDFM.cpp
  DoFManager.cpp
  OutputDataGPRS.cpp
  OutputDataVTK.cpp
  OutputDataPostprocessor.cpp
  VTKWriter.cpp
  MultiScaleDataMSRSB.cpp
  MultiScaleDataMech.cpp
  UnionFind.cpp
  DiscretizationFEM.cpp
)

SET_TARGET_PROPERTIES (
    preprocessor
    PROPERTIES LINKER_LANGUAGE CXX
)

set(preprocessor_libs ${preprocessor_libs} angem muparser mesh yaml-cpp)

if(METIS_FOUND)
  set(preprocessor_libs ${preprocessor_libs} ${METIS_LIBRARIES})
endif()

TARGET_LINK_LIBRARIES(preprocessor ${preprocessor_libs})
