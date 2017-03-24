/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
program SpinChainSZ
  use ESSEXCRS
  implicit none

  !chain length
  integer L

  !periodic boundary conditions?
  logical usePBC

  !physical parameters
  !spin-spin interaction
  double precision Jxy,Jz

  !number of spins up, 0 <= NUp <=L
  integer NUp
  
  !number of states
  integer(kind=LongIntType) Ns
  
  !the Hamilton operator
  type(SparseCRS) HMat

  integer(kind=LongIntType) :: il,jl,kl,bp,bp2,lbp,tbp,NE


  integer(kind=LongIntType), allocatable :: imask(:)


  !lookup tables and counter
  integer(kind=LongIntType), allocatable :: lutlead(:),luttrail(:),cnt(:),lutrev(:),revcnt(:)

  integer :: j

  double precision :: hp


  character(len=200) str,fn_out


  !!!!!!!!!!!!!!!!!!!!!!!!
  ! model parameters: Jxy,Jz, L
  ! Sz-symmetry: preserves number of spins up = Nup
  !
  ! vector length: <= 2^L
  ! matrix size: <= (2L+1)*2^L non-zero entries
  !

  call GetArg(1,str)
  read(str,*) Jxy
  call GetArg(2,str)
  read(str,*) Jz
  call GetArg(3,str)
  read(str,*) L
  call GetArg(4,str)
  read(str,*) NUp
  call GetArg(5,fn_out)
  
  ! periodic boundary conditions -> true
  usePBC=.true.

  if (mod(L,2) /= 0) stop 'abort -> need even L!'

  if ((NUp<0).or.(NUp>L)) stop 'abort -> need 0 <= NUp <= L'

  Ns=Binomial(L,NUp)
  
  write(6,*) 'Number of states: ',NS


  !create lookup tables and counter
  !+1 in lutlead because of loop below
  allocate(lutlead(2**(L/2+1)),luttrail(2**(L/2)),cnt(0:L),lutrev(2**(L/2)),revcnt(0:L))

  cnt=0
  do il=1,2**(L/2)
     j=bitcount(il-1)
     cnt(j)=cnt(j)+1
     luttrail(il)=cnt(j)
  end do

  !we create yet another lookup table. It's small, anyway.
  revcnt=0
  revcnt(0)=1
  do j=1,L
     revcnt(j)=revcnt(j-1)+cnt(j-1)
  end do
  do bp=0,2**(L/2)-1
     j=bitcount(bp)
     lutrev(revcnt(j))=bp
     revcnt(j)=revcnt(j)+1
  end do
  !and restore (or keep & copy)
  revcnt(0)=1
  do j=1,L
     revcnt(j)=revcnt(j-1)+cnt(j-1)
  end do

  lutlead(1)=1
  do il=1,2**(L/2)
     j=bitcount(il-1)
     j=Nup-j
     if ((j>=0).and.(j<=Nup)) then 
        lutlead(il+1)=lutlead(il)+cnt(j)
        kl=kl+cnt(j)
     else
        lutlead(il+1)=lutlead(il)
     end if
  end do

  kl=lutlead(2**(L/2)+1)-1

  write(6,*) 'kl= ',kl

  
  !pattern ---11---
  allocate(imask(L)) 
  imask(1)=3
  do j=2,L
     imask(j)=ishftc(imask(j-1),1,L)
  end do
 

  !create CRS matrix 

  !allocate matrix

  HMat%NR=Ns
  HMat%NC=Ns

  allocate(HMat%rptr(HMat%NR+1))
  
  !at most 2L+1 non-zero entries per row
  HMat%NE=(L+1)*Ns
  write(6,*) 'expected NE: ',HMAT%NE
  !note: final NE can be (much?) smaller than expected NE
  
  !likely to waste of memory here
  allocate(HMat%cind(HMat%NE))
  allocate(HMat%val(HMat%NE))
  
  !no symmetry
  HMat%Symmetric=.false.


  !construct matrix

  NE=0

  do il=1,2**(L/2)
     if (mod(il,2**10)==0) print*,il
     !jl runs from 1 ... Ns
     do jl=lutlead(il),lutlead(il+1)-1

        !**************        
        !we create all entries in row jl
        !from here ...

        !find bit pattern of spin state
        lbp=il-1   !index starts at 1
        j=NUp-bitcount(lbp)    
        tbp=lutrev(jl-lutlead(il)+revcnt(j))
        bp=ishft(lbp,L/2) + tbp

        !check program logic
        if (bitcount(bp) /= Nup) stop 'Error'

        !fill one row
        HMat%rptr(jl)=NE+1    

        !diagonal term
        
        hp=0.d0
        
        !spin-spin interaction
        !z-z
        if (Jz /= 0.d0) then
           do j=1,L
              if ((j==L).and.(.not.usePBC)) exit
              
              kl=ishftc(iand(bp,imask(j)),-(j-1),L)
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
           HMat%cind(NE)=jl
           HMat%val(NE)=hp
        end if
        
        !off-diagonal terms
        
        !spin-spin interaction
        !x-x + y-y
        
        if (Jxy /= 0.d0) then
           do j=1,L
              if ((j==L).and.(.not.usePBC)) exit
              
              kl=ishftc(iand(bp,imask(j)),-(j-1),L)
              if ((kl == 1).or.(kl == 2)) then
                 
                 bp2=ieor(bp,imask(j))
                 
                 !look up bitpattern
                 lbp = ishft(bp2,-L/2)
                 tbp = iand(bp2,2**(L/2)-1)
        
                 !index starts at 1
                 kl=lutlead(lbp+1) + luttrail(tbp+1) -1
                 
                 NE=NE+1
                 HMat%cind(NE)=kl    !index starts at 1
                 HMat%val(NE)=2.d0*Jxy
              end if
           end do
        end if
        
        !... to here
        !**************
          
     end do

  end do
  
 
  HMat%rptr(HMat%NR+1) = NE+1
  HMat%NE=NE

  write(6,*) 'NE = ',NE
  
  call WriteCRS(fn_out,HMat)



contains


  function Binomial(N,K)
    integer N,K
    integer(kind=LongIntType) Binomial
    integer i
    Binomial=1
    !don't change execution order
    do i=1,k
       Binomial=Binomial*(n-k+i)/i
    end do
  end function Binomial
  
  function bitcount(i)
    integer(kind=LongIntType) i,ii
    integer bitcount
    bitcount=0
    ii=i
    do while (bitcount < ii)
       bitcount=bitcount+1
       ii=iand(ii,ii-1)
    end do
  end function bitcount

end program SpinChainSZ
