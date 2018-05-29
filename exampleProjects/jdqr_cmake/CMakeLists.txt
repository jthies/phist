cmake_minimum_required(VERSION 3.5.0)

# project name
project(Djdqr CXX)

# find phist 
find_package(phist REQUIRED)
include_directories(${PHIST_INCLUDE_DIRS})
link_libraries(${PHIST_LIBRARIES})

set(CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED On)
add_executable(Djdqr Djdqr.cpp)
