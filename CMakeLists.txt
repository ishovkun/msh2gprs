cmake_minimum_required (VERSION 2.6)

# For YouCompleteMe
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CXX_COMPILER g++)

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/angem)
ADD_SUBDIRECTORY(src/angem)

# muparser
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/muparser)
ADD_SUBDIRECTORY(src/muparser)

# simdata
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/gprs-data)
ADD_SUBDIRECTORY(src/gprs-data)

# parser
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/parsers)
ADD_SUBDIRECTORY(src/parsers)


# targets
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17")
ADD_EXECUTABLE(msh2gprs ${CMAKE_SOURCE_DIR}/src/main.cpp)


# last flag for std::filesystem
# TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers -lstdc++fs -lasan)
TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers -lstdc++fs)
