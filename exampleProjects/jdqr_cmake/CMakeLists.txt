cmake_minimum_required(VERSION 2.8.10)

# project name
project(Djdqr CXX)

# find phist 
find_package(phist REQUIRED)
include_directories(${PHIST_INCLUDE_DIRS})
link_libraries(${PHIST_LIBRARIES})

# we also need MPI, but this should be handled by phist!
#find_package(MPI REQUIRED)
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
#link_libraries(${MPI_CXX_LIBRARIES})

set(CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED On)
add_executable(Djdqr Djdqr.cpp)
