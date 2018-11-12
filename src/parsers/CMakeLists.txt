# parser/CMakeLists.txt

ADD_LIBRARY(parsers
  # HEADERS
  json.hpp
  Parser.hpp
  GmshReader.hpp
  # IMPLEMENTATION
  Parser.cpp
  GmshReader.cpp
)

SET_TARGET_PROPERTIES (
    parsers
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(parsers PUBLIC
	${CMAKE_SOURCE_DIR}/src
)

TARGET_LINK_LIBRARIES(parsers gprs_data angem mesh)
