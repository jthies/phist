cmake_minimum_required(VERSION 2.8.10)

# project name
project(Djdqr CXX)

# find phist with pkg-config
include(FindPkgConfig)
pkg_check_modules(PHIST REQUIRED phist)
include_directories(${PHIST_INCLUDE_DIRS})
link_directories(${PHIST_LIBRARY_DIRS})
link_libraries(${PHIST_LIBRARIES})

# we also need MPI
find_package(MPI REQUIRED)
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
link_libraries(${MPI_CXX_LIBRARIES})

add_executable(Djdqr Djdqr.cpp)
