if(PHIST_HAVE_MKL)
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -DMKL_LP64")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -DMKL_LP64")
endif()

# -ffast kills high precision stuff!
set (CMAKE_C_FLAGS_RELEASE        "-O4 -fast -Mvect=levels:10")
set (CMAKE_CXX_FLAGS_RELEASE      "-O4 -fast -Mvect=levels:10")
set (CMAKE_Fortran_FLAGS_RELEASE  "-O4 -fast -Mvect=levels:10")

set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -g")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -g")
set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")

set (CMAKE_C_FLAGS_DEBUG       "-O0 -g")#-O0 -g -Mbounds -Mdwarf3")
set (CMAKE_CXX_FLAGS_DEBUG     "-O0")#-O0 -g -Mbounds -Mdwarf3")
set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")#-O0 -g -Mbounds -Mdwarf3")
