# preprocessor/CMakeLists.txt

ADD_LIBRARY(preprocessor
  Preprocessor.cpp
  parsers/YamlParser.cpp
  gmsh_interface/GmshInterface.cpp
  gmsh_interface/element_types.cpp
  utils.cpp
  CellPropertyManager.cpp
  EmbeddedFractureManager.cpp
  DiscreteFractureManager.cpp
  BoundaryConditionManager.cpp
  bc_expressions/location_functions.cpp
  bc_expressions/value_functions.cpp
  property_expressions/value_functions.cpp
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
  intersections/GridIntersectionSearcher.cpp
  intersections/UniformCartesianGrid.cpp
  intersections/Mapper.cpp
  discretization/mechanics/DiscretizationFEM.cpp
  discretization/mechanics/PolyhedralElementBase.cpp
  discretization/mechanics/PolyhedralElementMSRSB.cpp
  discretization/mechanics/JacobiPreconditioner.cpp
  discretization/mechanics/PolyhedralElementDirect.cpp
  discretization/mechanics/StandardFiniteElement.cpp
  discretization/mechanics/FEMFaceDoFManager.cpp
  discretization/mechanics/FeValuesTriangle.cpp
  discretization/mechanics/FeValuesWedge.cpp
  discretization/mechanics/FeValuesTetrahedron.cpp
  discretization/mechanics/FeValuesHexahedron.cpp
  discretization/mechanics/FeValuesQuadrangle.cpp
  discretization/mechanics/PFEM_integration/TributaryRegion3dBase.cpp
  discretization/mechanics/PFEM_integration/TributaryRegion3dFaces.cpp
  discretization/mechanics/PFEM_integration/TributaryRegion3dVertices.cpp
  discretization/mechanics/PFEM_integration/IntegrationRule3dFull.cpp
  discretization/mechanics/PFEM_integration/TributaryRegion2dFaces.cpp
  discretization/mechanics/PFEM_integration/TributaryRegion2dVertices.cpp
  discretization/mechanics/PFEM_integration/IntegrationRule3dAverage.cpp
  discretization/mechanics/PFEM_integration/IntegrationRule2dAverage.cpp
  discretization/mechanics/PFEM_integration/IntegrationRule2dPointwise.cpp
  discretization/mechanics/PFEM_integration/IntegrationRule2dFull.cpp
  discretization/mechanics/PFEM_integration/IntegrationRuleFractureAverage.cpp
  discretization/mechanics/PFEM_integration/IntegrationRuleFractureFull.cpp
  discretization/mechanics/PFEM_integration/IntegrationRule3dPointwise.cpp
  discretization/mechanics/PFEM_integration/FaceSorter.cpp
)

SET_TARGET_PROPERTIES (
    preprocessor
    PROPERTIES LINKER_LANGUAGE CXX
)

set(preprocessor_libs ${preprocessor_libs} angem muparser mesh yaml-cpp logger)

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
