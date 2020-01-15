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
  # element.cpp
  OutputDataGPRS.cpp
  OutputDataVTK.cpp
  OutputDataPostprocessor.cpp
  VTKWriter.cpp
  # renum.cpp
  # simdata.cpp
  # FlowData.cpp
  # transes.cpp
  # Well.cpp
  # MultiScaleDataMSRSB.cpp
  # MultiScaleDataMech.cpp
  # UnionFind.cpp
)

SET_TARGET_PROPERTIES (
    preprocessor
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(preprocessor PUBLIC
  ${CMAKE_SOURCE_DIR}/src
)

set(preprocessor_libs ${preprocessor_libs} angem muparser mesh yaml-cpp)

if(METIS_FOUND)
  set(preprocessor_libs ${preprocessor_libs} ${METIS_LIBRARIES})
endif()

TARGET_LINK_LIBRARIES(preprocessor ${preprocessor_libs})
