find_path(MKL_INCLUDE_DIR mkl.h $ENV{MKLROOT}/include)
set(MKL_LIBS $ENV{MKL_LIB})
include_directories (${MKL_INCLUDE_DIR} )
set(LIBS ${MKL_LIBS} dl)

#add_executable(ex1.bin ex1.cpp)
add_executable(ex2.bin ex2.cpp)
#target_link_libraries(ex1.bin ${LIBS} )
target_link_libraries(ex2.bin ${LIBS} )
