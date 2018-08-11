# gprs-data/CMakeLists.txt

ADD_LIBRARY(gprs_data
  # HEADERS
  element.hpp
  femout.hpp
  renum.hpp
  simdata.hpp
  transes.hpp
  utils.hpp
  SimdataConfig.hpp
  # IMPLEMENTATION
  element.cpp
  femout.cpp
  renum.cpp
  simdata.cpp
  transes.cpp
)

SET_TARGET_PROPERTIES (
    gprs_data
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_LINK_LIBRARIES(gprs_data angem muparser)
