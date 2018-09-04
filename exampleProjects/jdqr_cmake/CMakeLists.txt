cmake_minimum_required(VERSION 3.8.0)

# project name
project(Djdqr CXX)

# find phist 
find_package(phist REQUIRED)
include_directories(${PHIST_INCLUDE_DIRS})
link_libraries(${PHIST_LIBRARIES})

add_executable(Djdqr Djdqr.cpp)
