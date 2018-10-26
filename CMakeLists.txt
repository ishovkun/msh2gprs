cmake_minimum_required (VERSION 3.7)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/angem)
ADD_SUBDIRECTORY(src/angem)

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/mesh)
ADD_SUBDIRECTORY(src/mesh)

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
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

ADD_EXECUTABLE(msh2gprs ${CMAKE_SOURCE_DIR}/src/main.cpp)

TARGET_INCLUDE_DIRECTORIES(msh2gprs PRIVATE
	${CMAKE_SOURCE_DIR}/src
)

if (CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(msh2gprs PRIVATE "-Wall")
endif()

# memory leak debugging
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

# for profiling memory leaks
# TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers -lstdc++fs -lasan)

# last flag for std::filesystem
TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers -lstdc++fs)
