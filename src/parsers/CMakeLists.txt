# parser/CMakeLists.txt

ADD_LIBRARY(parsers
  # HEADERS
  json.hpp
  Parser.hpp
  # IMPLEMENTATION
  Parser.cpp
)

SET_TARGET_PROPERTIES (
    parsers
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_LINK_LIBRARIES(parsers gprs_data angem)
