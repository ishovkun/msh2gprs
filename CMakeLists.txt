cmake_minimum_required (VERSION 2.6)

# For YouCompleteMe
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/angem)
ADD_SUBDIRECTORY(src/angem)

# parser
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/parsers)
ADD_SUBDIRECTORY(src/parsers)

# simdata
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/gprs-data)
ADD_SUBDIRECTORY(src/gprs-data)



# targets
set(CMAKE_BUILD_TYPE Debug)

ADD_EXECUTABLE(msh2gprs ${CMAKE_SOURCE_DIR}/src/main.cpp)

TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers)
