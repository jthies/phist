/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
module ESSEXCRS
  implicit none

  
  !! kind parameters
  integer, parameter :: ShortIntType=4,LongIntType=8
  
  type SparseCRS
     !number of rows
     integer(kind=LongIntType) NR
     !number of columns
     integer(kind=LongIntType) NC
     !number of non-zero entries
     integer(kind=LongIntType) NE
     !row pointers
     integer(kind=LongIntType), allocatable :: rptr(:)
     !column indices
     integer(kind=LongIntType), allocatable :: cind(:)
     !matrix entries
     double precision, allocatable :: val(:)
     !symmetry flag
     logical Symmetric
  end type SparseCRS


contains
  
 subroutine WriteCRS(fname,CRS)
    character(*) fname
    type(SparseCRS) CRS
    ! TO CHANGE
    integer :: funit=10
    integer :: i

  

    open(funit,file=fname, access='STREAM', convert='LITTLE_ENDIAN', status='REPLACE')

    ! write header
    ! Endian indicator (here: Little Endian)
    write(funit) 0_ShortIntType
    ! version data format
    write(funit) 1_ShortIntType
    ! FORTRAN: Start counting with !
    write(funit) 0_ShortIntType
    ! symmetry information 
    if (CRS%Symmetric) then
       ! real symmetric
       write(funit) 2_ShortIntType
    else
       ! general
       write(funit) 1_ShortIntType
    end if
    ! data type (here: real double)
    write(funit) 6_ShortIntType
    !number of rows
    write(funit) CRS%NR
    !number of columns
    write(funit) CRS%NC
    !number of nonzero entries
    write(funit) CRS%NE
    !row pointers
    do i=1,CRS%NR+1
      write(funit) CRS%RPTR(i)-1
    end do
    !column indices
    do i=1,CRS%NE
      write(funit) CRS%CIND(i)-1
    end do
    !values
    write(funit) CRS%VAL(1:CRS%NE)
    close(funit)
    

  end subroutine WriteCRS
  
end module ESSEXCRS
