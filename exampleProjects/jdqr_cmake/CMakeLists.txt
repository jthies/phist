cmake_minimum_required(VERSION 3.8.0)

# project name
project(Djdqr CXX)

# find phist 
find_package(phist REQUIRED)

set(CMAKE_C_COMPILER ${PHIST_C_COMPILER})
set(CMAKE_CXX_COMPILER ${PHIST_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${PHIST_Fortran_COMPILER})

set(CMAKE_C_FLAGS "${PHIST_C_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${PHIST_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${PHIST_Fortran_FLAGS} ${CMAKE_CXX_FLAGS}")

include_directories(${PHIST_INCLUDE_DIRS})
link_directories(${PHIST_LIBRARY_DIR})

add_executable(Djdqr Djdqr.cpp)
target_link_libraries(Djdqr ${PHIST_LIBRARIES})
