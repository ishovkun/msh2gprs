# parser/CMakeLists.txt

# parser
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_DIR}/yaml)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/parsers/yaml)

ADD_LIBRARY(parsers
  # HEADERS
  json.hpp
  JsonParser.hpp
  YamlParser.hpp
  GmshReader.hpp
  # IMPLEMENTATION
  JsonParser.cpp
  YamlParser.cpp
  GmshReader.cpp
)

SET_TARGET_PROPERTIES (
    parsers
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(parsers PUBLIC
	${CMAKE_SOURCE_DIR}/src
	${CMAKE_SOURCE_DIR}/src/parsers/yaml/include
)

TARGET_LINK_LIBRARIES(parsers gprs_data angem mesh yaml-cpp)
