find_package(PythonInterp REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
file(GLOB INPUT_FILES ${PROJECT_SOURCE_DIR}/src/*/*.h*)

set(STATIC_HPP_FILES phist_types.hpp phist_exception.hpp)

set(HPP_FILES phist_kernels.hpp
phist_kernels_decl.hpp
phist_core.hpp
phist_core_decl.hpp
phist_jada.hpp
phist_jada_decl.hpp
phist_krylov.hpp
phist_krylov_decl.hpp
phist_precon.hpp
phist_precon_decl.hpp)

foreach (hpp_file ${STATIC_HPP_FILES})
  configure_file(${PROJECT_SOURCE_DIR}/cxx_bindings/${hpp_file} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()

set(ALL_HPP_FILES ${HPP_FILES} ${STATIC_HPP_FILES})

# put correct path into py-script and copy it
configure_file(${PROJECT_SOURCE_DIR}/cxx_bindings/make_cxx_bindings.py.in ${CMAKE_CURRENT_BINARY_DIR}/make_cxx_bindings.py)

# generate the C++ headers for all phist modules (subdirectories)
  add_custom_command(
      OUTPUT ${HPP_FILES}
      COMMAND  ${PYTHON_EXECUTABLE}
      ARGS ${CMAKE_CURRENT_BINARY_DIR}/make_cxx_bindings.py
      DEPENDS ${INPUT_FILES}
 )

add_custom_target(cxx_headers DEPENDS ${HPP_FILES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_cxx_bindings test_cxx_bindings.cpp)
target_link_libraries(test_cxx_bindings ${library_list})

add_dependencies(test_cxx_bindings cxx_headers)

enable_testing()
add_test(NAME cxx_bindings COMMAND ${TEST_ENV} ./test_cxx_bindings)
set_tests_properties(cxx_bindings PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*PASSED[. \n]*")
set_tests_properties(cxx_bindings PROPERTIES FAIL_REGULAR_EXPRESSION ". \n]*FAILED[. \n]*")
set_tests_properties(cxx_bindings PROPERTIES TIMEOUT 30)

# install all generated headers
foreach (hpp_file ${HPP_FILES})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${hpp_file} DESTINATION include)
endforeach()

foreach (hpp_file ${STATIC_HPP_FILES})
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${hpp_file} DESTINATION include)
endforeach()

