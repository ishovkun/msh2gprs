# gprs-data/CMakeLists.txt

ADD_LIBRARY(gprs_data
  element.cpp
  femout.cpp
  renum.cpp
  simdata.cpp
  FlowData.cpp
  transes.cpp
  VTKWriter.cpp
  Well.cpp
  MultiScaleDataMSRSB.cpp
)

SET_TARGET_PROPERTIES (
    gprs_data
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(gprs_data PUBLIC
	${CMAKE_SOURCE_DIR}/src
)

TARGET_LINK_LIBRARIES(gprs_data angem muparser mesh metis)
