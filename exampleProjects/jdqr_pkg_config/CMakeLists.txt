cmake_minimum_required(VERSION 2.8.10)

# project name
project(Djdqr CXX)

# find phist with pkg-config
find_package(PkgConfig)
pkg_check_modules(PHIST REQUIRED phist)
# the FindPkgConfig module creates a list rather than a string,
# for the CFLAGS, I think the must be a more elegant way than this:
string(REPLACE ";" " " PHIST_CFLAGS_STRING "${PHIST_CFLAGS_OTHER}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PHIST_CFLAGS_STRING}") 
include_directories(${PHIST_INCLUDE_DIRS})
link_directories(${PHIST_LIBRARY_DIRS})
link_libraries(${PHIST_LIBRARIES})

if (NOT PHIST_LIBRARIES MATCHES "phist_kernels_eigen")
  # we also need MPI
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
  link_libraries(${MPI_CXX_LIBRARIES})
endif()

add_executable(Djdqr Djdqr.cpp)
