# msh2gprs_root/CMakeLists.txt

cmake_minimum_required (VERSION 3.8)

project(ADGPRS_PREPROCESSOR)

set(CMAKE_CXX_STANDARD 17)

# for ycmd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

# boost is optional for 256-bit integers (angem point set)
# but it works muuuch faster with it
find_package(Boost 1.65.1)
if(Boost_FOUND)
  message("Boost found. Using boost libraries for hashing.")
  message("To turn off edit the CMakeLists file.")
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
  add_definitions(-DUSE_BOOST)
else()
  # 256 bit integer lib (for mesh)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/external/uint256)
  ADD_SUBDIRECTORY(external/uint256)
ENDIF()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
find_package(METIS)
if(METIS_FOUND)
  set(GLOBAL_INCLUDE_DIRS ${GLOBAL_INCLUDE_DIRS} ${METIS_INCLUDE_DIRS})
  message("-- METIS found")
  add_definitions(-DWITH_METIS)
else()
  message("-- METIS not found")
endif()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/external)

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/external/angem)
ADD_SUBDIRECTORY(external/angem)
INCLUDE_DIRECTORIES(${angem_include_dirs})

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/mesh)
ADD_SUBDIRECTORY(src/mesh)
# muparser
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/external/muparser)
ADD_SUBDIRECTORY(external/muparser)
# yaml parser
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/external/yaml)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/external/yaml)

# preprocessor
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/preprocessor)
ADD_SUBDIRECTORY(src/preprocessor)

# targets
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

ADD_EXECUTABLE(msh2gprs ${CMAKE_SOURCE_DIR}/src/main.cpp)

TARGET_INCLUDE_DIRECTORIES(msh2gprs PRIVATE
  ${CMAKE_SOURCE_DIR}/external
  ${CMAKE_SOURCE_DIR}/src
)

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
endif()

# # memory leak debugging
# # set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# # set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# for profiling memory leaks
# TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers mesh uint256 -lstdc++fs -lasan)

# last flag for std::filesystem
TARGET_LINK_LIBRARIES(msh2gprs preprocessor yaml-cpp mesh -lstdc++fs)

enable_testing()
add_subdirectory(tests)
