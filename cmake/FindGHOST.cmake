# Copyright (c) 2012 Marwan Abdellah <marwan.abdellah@epfl.ch>
#                    Daniel Nachbaur <daniel.nachbaur@epfl.ch>
#               2013 Stefan.Eilemann@epfl.ch

# Use pkg-config to fetch the contents of the .pc file
# After that, use the directories refer to the libraries and
# also the headers

find_package(PkgConfig)

if(GHOST_ROOT)
  set(ENV{PKG_CONFIG_PATH} "${GHOST_ROOT}/lib/pkgconfig")
else()
  foreach(PREFIX ${CMAKE_PREFIX_PATH})
    set(PKG_CONFIG_PATH "${PKG_CONFIG_PATH}:${PREFIX}/lib/pkgconfig")
  endforeach()
  set(ENV{PKG_CONFIG_PATH} "${PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
endif()

if(ghost_FIND_REQUIRED)
  set(_ghost_OPTS "REQUIRED")
elseif(ghost_FIND_QUIETLY)
  set(_ghost_OPTS "QUIET")
else()
  set(_ghost_output 1)
endif()

if(ghost_FIND_VERSION)
  if(ghost_FIND_VERSION_EXACT)
    pkg_check_modules(GHOST ${_ghost_OPTS} ghost=${ghost_FIND_VERSION})
  else()
    pkg_check_modules(GHOST ${_ghost_OPTS} ghost>=${ghost_FIND_VERSION})
  endif()
else()
  pkg_check_modules(GHOST ${_ghost_OPTS} ghost)
endif()

if(GHOST_FOUND)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GHOST DEFAULT_MSG GHOST_LIBRARIES
    GHOST_INCLUDE_DIRS GHOST_LIBRARY_DIRS)

  if(NOT ${GHOST_VERSION} VERSION_LESS 1.7.0)
    set(GHOST_GL_FOUND 1)
  endif()

  if(_ghost_output)
    message(STATUS
      "Found ghost ${GHOST_VERSION} in ${GHOST_INCLUDE_DIRS}:${GHOST_LIBRARIES}")
  endif()
endif()
