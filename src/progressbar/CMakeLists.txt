# progressbar/CMakeLists.txt

ADD_LIBRARY(progressbar
  src/ProgressBar.cpp)

SET_TARGET_PROPERTIES (
    progressbar
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(progressbar PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include/progressbar
  )

set(progressbar_include_dirs ${progressbar_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR}/include  PARENT_SCOPE)
