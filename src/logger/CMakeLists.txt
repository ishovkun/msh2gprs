# progressbar/CMakeLists.txt

ADD_LIBRARY(logger
  src/Color.cpp
  src/ProgressBar.cpp
  src/Logger.cpp
  )

SET_TARGET_PROPERTIES (
  logger
  PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(logger PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include/logger
  )

set(logger_include_dirs ${logger_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR}/include  PARENT_SCOPE)
