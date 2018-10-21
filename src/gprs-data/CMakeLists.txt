# gprs-data/CMakeLists.txt

ADD_LIBRARY(gprs_data
  # HEADERS
  element.hpp
  femout.hpp
  renum.hpp
  transes.hpp
  simdata.hpp
  utils.hpp
  SimdataConfig.hpp
  FlowData.hpp
  VTKWriter.hpp
  # IMPLEMENTATION
  element.cpp
  femout.cpp
  renum.cpp
  simdata.cpp
  FlowData.cpp
  transes.cpp
  VTKWriter.cpp
)

SET_TARGET_PROPERTIES (
    gprs_data
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_LINK_LIBRARIES(gprs_data angem muparser)
