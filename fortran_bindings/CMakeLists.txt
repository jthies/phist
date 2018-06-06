#cmake_minimum_required(VERSION 3.8.0)
project(phist_fort Fortran)

enable_testing()

set(PHIST_SRC ${PROJECT_SOURCE_DIR}/../src)
set(INC_FLAGS "-I${PHIST_SRC}/kernels/"
              "-I${PHIST_SRC}/kernels/${PHIST_KERNEL_LIB}/"
              "-I${PHIST_SRC}/tools/"
              "-I${PHIST_SRC}/precon/"
              "-I${PHIST_SRC}/jada/"
              "-I${PHIST_SRC}/krylov/"
              "-I${PHIST_SRC}/core/"
              "-I${PROJECT_SOURCE_DIR}/../drivers/"
              "-I${CMAKE_CURRENT_BINARY_DIR}/../"
              )

set(F_EXT "F90" CACHE STRING "file extension to be added to generated Fortran modules (.F90 by default)")
set(MAX_LINE_LEN "120" CACHE STRING "break lines in Fortran output after a certain number of characters")

# modules independent of data type
set(BASENAMES_TI "kernels"
              "jada"
              "tools"
              )

# typed modules
set(BASENAMES "tools"
              "kernels"
              "core"
              "krylov"
              "jada"
              "precon"
              )

set(TEST_DRIVERS "kernels"
                 "core"
                 "jada"
              )

# dynamically generated Fortran source files
set(FFILES "")
# static source file we need for index types
set(FSTATIC "${CMAKE_BINARY_DIR}/phist_types.F90"
            "${PROJECT_SOURCE_DIR}/phist_testing.F90")

set(CMAKE_Fortran_FLAGS "-I${PROJECT_SOURCE_DIR} ${CMAKE_Fortran_FLAGS}")

############################
# type-independent modules #
############################

configure_file(${PROJECT_SOURCE_DIR}/phist_types.F90 ${CMAKE_BINARY_DIR}/phist_types.F90 COPYONLY)

foreach (base ${BASENAMES_TI})
  
  set(TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/${base})

  list(APPEND HFILES "${TARGET_DIR}/${base}.h")
  list(APPEND FFILES "phist_${base}.${F_EXT}")
  list(APPEND OFILES "phist_${base}.o")

# copy input file and include type-specific headers and CMake options
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${base}.h.in" "${TARGET_DIR}/${base}-1.CFWRAPPER_SKIP.h")

# preprocess C header file
  add_custom_command(
      OUTPUT ${TARGET_DIR}/${base}-2.CFWRAPPER_SKIP.h
      COMMAND  gcc
      ARGS ${INC_FLAGS} ${C_DEFINES} ${C_FLAGS} -E ${base}/${base}-1.CFWRAPPER_SKIP.h > ${base}/${base}-2.CFWRAPPER_SKIP.h
      DEPENDS ${TARGET_DIR}/${base}-1.CFWRAPPER_SKIP.h
 )

# strip of preprocessor comments
    add_custom_command(
        OUTPUT ${TARGET_DIR}/${base}.h
        COMMAND  sed
        ARGS  "s/\"^#.*\"//g" ${TARGET_DIR}/${base}-2.CFWRAPPER_SKIP.h > ${TARGET_DIR}/${base}.h
        DEPENDS ${TARGET_DIR}/${base}-2.CFWRAPPER_SKIP.h
   )
endforeach()
#
# always generate all interfaces for now, because we otherwise
# have to adjust the cfwrapper file
set(TG_Types S D C Z)

foreach (base ${BASENAMES})
  foreach (tc ${TG_Types})
    string(TOLOWER ${tc} tcl)
    set(TG_HFILE phist_gen_${tcl}.h)

    set(tcl_base ${base}_${tcl})
    set(TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/${tcl_base})

    list(APPEND HFILES "${TARGET_DIR}/${tcl_base}.h")
  list(APPEND FFILES "phist_${tcl_base}.${F_EXT}")
  list(APPEND OFILES "phist_${tcl_base}.o")

# copy input file and include type-specific headers and CMake options
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${base}_X.h.in" "${TARGET_DIR}/${tcl_base}-1.CFWRAPPER_SKIP.h")

# preprocess C header file
    add_custom_command(
        OUTPUT ${TARGET_DIR}/${tcl_base}-2.CFWRAPPER_SKIP.h
        COMMAND  gcc
        ARGS ${INC_FLAGS} ${C_DEFINES} ${C_FLAGS} -E ${tcl_base}/${tcl_base}-1.CFWRAPPER_SKIP.h > ${tcl_base}/${tcl_base}-2.CFWRAPPER_SKIP.h
        DEPENDS ${TARGET_DIR}/${tcl_base}-1.CFWRAPPER_SKIP.h
   )
# strip of preprocessor comments
    add_custom_command(
        OUTPUT ${TARGET_DIR}/${tcl_base}.h
        COMMAND  sed
        ARGS  "s/\"^#.*\"//g" ${TARGET_DIR}/${tcl_base}-2.CFWRAPPER_SKIP.h > ${TARGET_DIR}/${tcl_base}.h
        DEPENDS ${TARGET_DIR}/${tcl_base}-2.CFWRAPPER_SKIP.h
   )

  endforeach()
endforeach()

set(DICT phist_dict)
configure_file("${PROJECT_SOURCE_DIR}/phist_dict.py.in" "${CMAKE_CURRENT_BINARY_DIR}/phist_dict.py")
configure_file("${PROJECT_SOURCE_DIR}/cfwrapper.py.in" "${CMAKE_CURRENT_BINARY_DIR}/cfwrapper.py")

add_custom_command(
        OUTPUT ${FFILES}
        COMMAND  ${CMAKE_CURRENT_BINARY_DIR}/cfwrapper.py
        DEPENDS ${HFILES}
   )

# the following phony targets are meant for debugging/development only

# generate the preprocessed header files
add_custom_target(pp DEPENDS ${HFILES} COMMAND ${CMAKE_COMMAND} -E remove */*.CFWRAPPER_SKIP.h)

# generate the preprocessed header files, keep temporary files
add_custom_target(pp_keep DEPENDS ${HFILES})

# run the cfwrapper.py script to create all fortran source files
add_custom_target(cf DEPENDS pp_keep COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cfwrapper.py)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_library(phist_fort ${FFILES} ${FSTATIC})


foreach(TESTNAME ${TEST_DRIVERS})
        add_executable(test_fortran_${TESTNAME} test/${TESTNAME}.F90)

        target_link_libraries(test_fortran_${TESTNAME} phist_fort ${library_list})

        add_test(NAME fortran_${TESTNAME} COMMAND ./test_fortran_${TESTNAME})
        set_tests_properties(fortran_${TESTNAME} PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*ALL TESTS PASSED[. \n]*")
        set_tests_properties(fortran_${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "[. \n]*FAILED[. \n]*")
        set_tests_properties(fortran_${TESTNAME} PROPERTIES TIMEOUT 500)
endforeach()

# top-level target for PHIST
add_custom_target(fortran_modules DEPENDS cf phist_fort)

