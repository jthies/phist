cmake_minimum_required(VERSION 2.8.10)

# project name and language. This is actually a C++ project
# but it loads phist. If phist was built using the builtin kernels,
# they depend on the Fortran variant of MPI, and if we don't claim
# here that this is a Fortran project, too, find_package(MPI) below will
# only give us the MPI_CXX variables. CMake projects don't have this issue.
#
# Obviously it would be nicer to handle the dependencies via pkg-config, but
# since there is no pkg-config file for MPI I can't add it as a dependency.
project(Djdqr CXX Fortran)

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
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
  # the Fortran_LIBRARIES are used by the builtin kernels only
  link_libraries(${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif()

add_executable(Djdqr Djdqr.cpp)
