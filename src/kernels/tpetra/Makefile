
-include ${ESSEX_BASE}/make.inc/Makefile.${ESSEX_PLAT}

CXXFLAGS+=-I${TRILINOS_BASE}/include

OBJ=kernels.o

LIBNAME=libessex_kernels_tpetra_${ESSEX_PLAT}.a
HEADER_FILES=../essex_kernels.h ../essex_kernels_decl.h ../gen_*.h

include ${ESSEX_BASE}/make.inc/rules
