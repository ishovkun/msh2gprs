# - Find Gmsh library
# Find the Gmsh includes and library
# This module defines
#  Gmsh_INCLUDE_DIR: where to find gmsh.h
#  Gmsh_LIBRARY: the Gmsh library
#  Gmsh_FOUND: if false, do not try to use Gmsh
#  Gmsh_VERSION: The found version of Gmsh

if (DEFINED GMSH_INSTALL_PATH)
  find_path(Gmsh_INCLUDE_DIR NAMES Gmsh.h gmsh.h
    PATHS
    ${GMSH_INSTALL_PATH}/include
    PATH_SUFFIXES gmsh
    NO_DEFAULT_PATH
    )
else()
  find_path(Gmsh_INCLUDE_DIR NAMES gmsh.h
  PATHS
  /usr/include/
  /usr/local/include/
  PATH_SUFFIXES
  gmsh
)
endif()

if (DEFINED GMSH_INSTALL_PATH)
  message ("searching ${GMSH_INSTALL_PATH}/lib")
  find_library(Gmsh_LIBRARY NAMES gmsh libgmsh libGmsh.so
    PATHS
    ${GMSH_INSTALL_PATH}/lib
    NO_DEFAULT_PATH
    )
else()
  find_library(Gmsh_LIBRARY NAMES gmsh libgmsh
  PATHS
  /usr/lib
  /usr/local/lib
  /usr/lib64
  /usr/local/lib64
  PATH_SUFFIXES
  gmsh
)
endif()

# get version
if(Gmsh_INCLUDE_DIR)
  set(Gmsh_INCLUDE_DIRS ${Gmsh_INCLUDE_DIR})

  if(EXISTS ${Gmsh_INCLUDE_DIR}/gmsh.h)
    set(HeaderFile gmsh.h)
  elseif(EXISTS ${Gmsh_INCLUDE_DIR}/Gmsh.h)
    set(HeaderFile Gmsh.h)
  endif()

  # set(VersionFile gmsh.h)
  if(EXISTS ${Gmsh_INCLUDE_DIR}/${HeaderFile})
    file(READ ${Gmsh_INCLUDE_DIR}/${HeaderFile} GMSH_VERSION_FILE)
    string(REGEX MATCH "\#define GMSH_API_VERSION *\"([0-9,.]*).*\"" GMSH_VERSION_STRING ${GMSH_VERSION_FILE})
    set(Gmsh_VERSION ${CMAKE_MATCH_1} CACHE INTERNAL "Gmsh Version")
  else()
    message(SEND_ERROR "Could not find " ${VersionFile} " in " ${Gmsh_INCLUDE_DIR}
                       ". Check path or reconfigure Gmsh with -DENABLE_PRIVATE_API=ON")
  endif()

  if(Gmsh_LIBRARY)
    set(Gmsh_LIBRARIES ${Gmsh_LIBRARY})
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set Gmsh_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  FOUND_VAR Gmsh_FOUND
  REQUIRED_VARS Gmsh_LIBRARY Gmsh_INCLUDE_DIR
  VERSION_VAR Gmsh_VERSION)

if(Gmsh_FOUND)
  set(GMSH_LIBRARIES ${Gmsh_LIBRARY})
  set(GMSH_INCLUDE_DIRS ${Gmsh_INCLUDE_DIR})
endif()

mark_as_advanced(
  Gmsh_INCLUDE_DIR
  Gmsh_LIBRARY
)
