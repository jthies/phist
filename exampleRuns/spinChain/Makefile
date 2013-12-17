FC=gfortran
FFLAGS=-O2 -g
src1=SpinChainCRS_v191113.f90
src2=SpinChainSZ_v101213.f90

default: SpinChain SpinChainSZ

essexcrs.mod: essexcrs.f90
	${FC} ${FFLAGS} -c $<

SpinChain: ${src1} essexcrs.mod essexcrs.o
	${FC} ${FFLAGS} -o $@ $< essexcrs.o

SpinChainSZ: ${src2} essexcrs.mod essexcrs.o
	${FC} ${FFLAGS} -o $@ $< essexcrs.o
	
