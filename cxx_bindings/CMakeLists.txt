find_package(PythonInterp REQUIRED)

file(GLOB INPUT_FILES ${PROJECT_SOURCE_DIR}/src/*/*.h*)

set(HPP_FILES phist_kernels.hpp
phist_kernels_impl.hpp
phist_core.hpp
phist_core_impl.hpp
phist_jada.hpp
phist_jada_impl.hpp
phist_krylov.hpp
phist_krylov_impl.hpp
phist_precon.hpp
phist_precon_impl.hpp)

# generate the C++ headers for all phist modules (subdirectories)
  add_custom_command(
      OUTPUT ${HPP_FILES}
      COMMAND  ${PYTHON_EXECUTABLE}
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/make_cxx_bindings.py
      DEPENDS ${INPUT_FILES}
 )


add_executable(test_cxx_bindings test_cxx_bindings.cpp)
#target_link_libraries(${PHIST_LIBRARIES})
