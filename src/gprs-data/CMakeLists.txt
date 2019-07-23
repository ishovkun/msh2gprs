# gprs-data/CMakeLists.txt

ADD_LIBRARY(gprs_data
  element.cpp
  femout.cpp
  OutputDataVTK.cpp
  renum.cpp
  simdata.cpp
  FlowData.cpp
  transes.cpp
  VTKWriter.cpp
  Well.cpp
  MultiScaleDataMSRSB.cpp
  MultiScaleDataMech.cpp
  UnionFind.cpp
)

SET_TARGET_PROPERTIES (
    gprs_data
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(gprs_data PUBLIC
  ${CMAKE_SOURCE_DIR}/src
)

set(grps_data_libs ${grps_data_libs} angem muparser mesh)

if(METIS_FOUND)
  set(grps_data_libs ${grps_data_libs} ${METIS_LIBRARIES})
endif()

TARGET_LINK_LIBRARIES(gprs_data ${grps_data_libs})
