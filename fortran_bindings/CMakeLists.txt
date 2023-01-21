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

if (PHIST_KERNEL_LIB_GHOST)
  list(GET GHOST_INCLUDE_DIRS 0 GHOST_INCLUDE_DIR)
  set(INC_FLAGS "-I${GHOST_INCLUDE_DIR}" ${INC_FLAGS})
elseif (PHIST_KERNEL_LIB_EPETRA)
  set(INC_FLAGS "-I${Epetra_DIR}/../../../include" ${INC_FLAGS})
elseif (PHIST_KERNEL_LIB_TPETRA)
  set(INC_FLAGS "-I${Tpetra_DIR}/../../../include" ${INC_FLAGS})
endif()

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
set(FSTATIC "${PROJECT_SOURCE_DIR}/phist_types.F90"
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
      COMMAND  ${CMAKE_C_COMPILER}
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
        COMMAND  ${CMAKE_C_COMPILER}
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
#add_custom_target(pp DEPENDS ${HFILES} COMMAND ${CMAKE_COMMAND} -E remove */*.CFWRAPPER_SKIP.h)

# generate the preprocessed header files, keep temporary files
#add_custom_target(pp_keep DEPENDS ${HFILES})

# run the cfwrapper.py script to create all fortran source files
#add_custom_target(cf DEPENDS pp_keep COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cfwrapper.py)

add_library(phist_fort ${FFILES} ${FSTATIC})
#set_target_properties(phist_fort PROPERTIES DEPENDS cf)
set_target_properties(phist_fort PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(phist_fort PUBLIC ${MPI_Fortran_MODULE_DIR})

# the type-specific modules depend on the general ones, e.g. phist_jada_s.F90 needs phist_jada.mod to compile.
# I think CMake should be able to figure out module dependencies, but maybe it's confused by our code generation.
foreach (base ${BASENAMES})
  foreach (tc ${TG_Types})
    set_source_files_properties(phist_${tcl_base}.${F_EXT} PROPERTIES DEPENDS phist_${base}.o)
  endforeach()
endforeach()

foreach(TESTNAME ${TEST_DRIVERS})
        add_executable(test_fortran_${TESTNAME} test/${TESTNAME}.F90)

        target_link_libraries(test_fortran_${TESTNAME} phist_fort ${library_list} ${MPI_Fortran_LIBRARIES})
        target_include_directories(test_fortran_${TESTNAME} PUBLIC ${MPI_Fortran_MODULE_DIR})

        add_test(NAME fortran_${TESTNAME} COMMAND ${TEST_ENV} ./test_fortran_${TESTNAME})
        set_tests_properties(fortran_${TESTNAME} PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*ALL TESTS PASSED[. \n]*")
        set_tests_properties(fortran_${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "[. \n]*FAILED[. \n]*")
        set_tests_properties(fortran_${TESTNAME} PROPERTIES TIMEOUT 500)
endforeach()

# makes sure headers and module files are found when building the f-files and tests, resp.
include_directories(${CMAKE_CURRENT_BINARY_DIR})


# top-level target for PHIST
add_custom_target(fortran_modules DEPENDS phist_fort)

# install all generated Fortran files. Note that we install the actual source files,
# not the modules, because I don't know how to get their names here. So when compiling
# a program using the modules, they will be regenerated.
foreach (f_file ${FFILES})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${f_file} DESTINATION include/)
endforeach()

foreach (f_file ${FSTATIC})
    install(FILES ${f_file} DESTINATION include/)
endforeach()

