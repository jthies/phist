program SpinChainCRS
  use ESSEXCRS
  implicit none

  !chain length
  integer L

  !periodic boundary conditions?
  logical usePBC

  !physical parameters
  !spin-spin interaction
  double precision Jx,Jy,Jz
  !magnetic field
  double precision Bx,Bz
  
  !number of states
  integer(kind=LongIntType) Ns
  
  !the Hamilton operator
  type(SparseCRS) HMat

  integer(kind=LongIntType) :: il,kl,NE


  integer(kind=LongIntType), allocatable :: imask(:)

  integer :: j

  double precision :: hp


  character(len=200) str,fn_out


  !!!!!!!!!!!!!!!!!!!!!!!!
  ! model parameters: Jx,Jy,Jz, Bx,By, L
  !
  ! vector length: 2^L
  ! matrix size: <= (2L+1)*2^L non-zero entries
  !

  call GetArg(1,str)
  read(str,*) Jx
  call GetArg(2,str)
  read(str,*) Jy
  call GetArg(3,str)
  read(str,*) Jz
  call GetArg(4,str)
  read(str,*) Bx
  call GetArg(5,str)
  read(str,*) Bz
  call GetArg(6,str)
  read(str,*) L
  call GetArg(7,fn_out)
  
  ! periodic boundary conditions -> true
  usePBC=.true.

  !
  !!!!!!!!!!!!!!!!!!!!!!!!


  allocate(imask(L)) 
  imask(1)=3
  do j=2,L
     imask(j)=ishftc(imask(j-1),1,L)
  end do
 
  Ns=2**L
  
  HMat%NR=Ns
  HMat%NC=Ns

  allocate(HMat%rptr(HMat%NR+1))

  !at most 2L+1 non-zero entries per row
  HMat%NE=(2*L+1)*Ns   !tentatively

  allocate(HMat%cind(HMat%NE))
  allocate(HMat%val(HMat%NE))

  !no symmetry used
  HMat%Symmetric=.false.


  NE=0
  
  !loop over all states,
  !represented as bit patterns with length L
  do il=0,Ns-1

     !fill one row
     HMat%rptr(il+1)=NE+1      !index starts at 1

     !diagonal term

     hp=0.d0

     !magnetic field
     if (Bz /= 0.d0) then
        do j=1,L
           if (btest(il,j-1)) then
              hp=hp+Bz
           else
              hp=hp-Bz
           end if
        end do
     end if
     
     !spin-spin interaction
     !z-z
     if (Jz /= 0.d0) then
        do j=1,L
           if ((j==L).and.(.not.usePBC)) exit
           
           kl=ishftc(iand(il,imask(j)),-(j-1),L)
           select case (kl)
           case (0)
              hp=hp+Jz
           case (1)
              hp=hp-Jz
           case (2)
              hp=hp-Jz
           case (3)
              hp=hp+Jz
           case default
              stop 'error 1'
           end select

        end do
     end if
     
     if (hp /= 0.d0) then
        NE=NE+1
        HMat%cind(NE)=il+1    !index starts at 1
        HMat%val(NE)=hp
     end if

     
     !off-diagonal terms

     !magnetic field
     if (Bx /= 0.d0) then
        do j=1,L
           kl=ieor(il,2**(j-1))
           NE=NE+1
           HMat%cind(NE)=kl+1   !index starts at 1
           HMat%val(NE)=Bx
        end do
     end if

 

     !spin-spin interaction
     !x-x
     !y-y
     
     if ((Jx /= 0.d0).or.(Jy /= 0.d0)) then
        do j=1,L
           if ((j==L).and.(.not.usePBC)) exit

           !unnecessary overhead (see previous loop)
           kl=ishftc(iand(il,imask(j)),-(j-1),L)
           select case (kl)
           case (0)
              hp=Jx-Jy
           case (1)
              hp=Jx+Jy
           case (2)
              hp=Jx+Jy
           case (3)
              hp=Jx-Jy
           case default
              stop 'error 2'
           end select  

           if (hp /= 0.d0) then
              kl=ieor(il,imask(j))
              NE=NE+1
              HMat%cind(NE)=kl+1    !index starts at 1
              HMat%val(NE)=hp
           end if
        end do
     end if

     
  end do

  
  HMat%rptr(HMat%NR+1) = NE+1
  HMat%NE=NE
  
  write(6,*) 'square matrix dimension: ',NS
  write(6,*) 'non-zero matrix entries: ',NE
    
  call WriteCRS(fn_out,HMat)



end program SpinChainCRS
