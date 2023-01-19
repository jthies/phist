#
# phist integration tests run the ./drivers with certain test problems and parameters to check that
# the iterative solvers converge as expected.
#

set(MPIRUN_COMMAND ${MPIEXeC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_POSTFLAGS} ${TEST_ENV})

# Jacobi-Davidson solver
add_test(NAME Dsubspacejada COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dsubspacejada  spinSZ14 I 1 10 SR 1e-8 28 4 16 32 4 8 0.0 0 1 1)


# iterative solvers: BiCGstab and IDR(s)
add_test(NAME Dbicgstab COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedidrs matpde64 0 288 1e-6)

add_test(NAME Didr1 COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedidrs matpde64 1 288 1e-6)

add_test(NAME Didr2 COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedidrs matpde64 2 262 1e-6)

add_test(NAME Didr4 COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedidrs matpde64 4 234 1e-6)

# GMRES
add_test(NAME Dgmres COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedgmres matpde64 150 1e-6 1 100)

add_test(NAME Dgmres_restarted COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedgmres DjadaTestMat512.mm 150 1e-6 1 30)

if (PHIST_HAVE_IFPACK OR PHIST_HAVE_IFPACK2)
    file(COPY ${CMAKE_SOURCE_DIR}/exampleRuns/precon/ifpack_options.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    # note that we're rather generous here because ifpack2 does (>60 sequentially) does worse than ifpack (42 iterations)
    # and the test should also work with more threads
    add_test(NAME Dgmres_ilu  COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_Dblockedgmres DjadaTestMat512.mm 100 1e-10 1 10 ifpack ifpack_options.xml)
endif()

if (PHIST_KERNEL_LIB_BUILTIN OR PHIST_KERNEL_LIB_EPETRA)

    # note that we're rather generous here because the way the KACZ/CARP kernel is implemented has some influence on the convergence speed.
    # and the test should also work with more MPI processes and different variants.
    file(COPY ${CMAKE_SOURCE_DIR}/test/matrices/graphene/128x64/A.mm DESTINATION {CMAKE_CURRENT_BINARY_DIR}/graphene128x64.mm)
    add_test(NAME Dcarp_cg_graphene COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_carp_cg graphene128 none 1 1 1e-12 350 1.0)

    file(COPY ${PROJECT_SOURCE_DIR}/exampleRuns/imaginary_shifts.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    # CARP-CG with complex shifts
    add_test(NAME Dcarp_cg_spin_complex_shifts COMMAND ${MPIRUN_COMMAND}
        ./drivers/phist_carp_cg spinSZ14 imaginary_zshifts.txt 1 1 1e-6 350 1.2')
endif()
