
-include ${ESSEX_BASE}/make.inc/Makefile.${ESSEX_PLAT}

CXXFLAGS+=-I${TRILINOS_BASE}/include

OBJ=kernels.o

LIBNAME=libessex_kernels_epetra_${ESSEX_PLAT}.a
HEADER_FILES=../kernels.h ../kernels_decl.h ../gen_*.h

${LIBNAME}: ${OBJ}
	${AR} ${ARFLAGS} $@ $<

include ${ESSEX_BASE}/make.inc/rules
