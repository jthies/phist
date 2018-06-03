module phist_testing

use mpi

public :: tests_init, tests_finalize, parallel_expect, parallel_assert

private 

character(len=128) :: current_test_file
logical :: all_tests_passed
integer :: comm, ierr


contains

subroutine tests_init(filename)
use mpi
implicit none
character(len=*) :: filename

comm=MPI_COMM_WORLD
current_test_file=filename
all_tests_passed=.true.

end subroutine tests_init

subroutine tests_finalize(printer)
implicit none
logical :: printer

if (all_tests_passed .and. printer) then
  write(*,*) 'ALL TESTS PASSED'
end if

end subroutine tests_finalize

subroutine parallel_expect(result,line,printer)
use mpi
implicit none
logical result, global_result, printer
integer :: line

call MPI_Allreduce(result,global_result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
all_tests_passed = all_tests_passed .and. global_result
if ((.not. global_result) .and. printer) then
        write(*,*) 'TEST FILE: ',current_test_file
        write(*,*) 'ASSERTION FAILED IN LINE ',line
end if
end subroutine parallel_expect

subroutine parallel_assert(result,line,printer)
use mpi
implicit none
logical result, printer
integer :: line

call parallel_expect(result,line,printer)
if (.not. result) then
  call MPI_Abort(comm,1,ierr)
end if
end subroutine parallel_assert

end module phist_testing
