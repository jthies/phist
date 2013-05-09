# this makefile requires that you set the environment variable
# "ESSEX_KERNEL_LIB" to either
#
# - "epetra"
# - "tpetra"
# - "ghost"
#
# and builds the according library called kernels.a
#

default: install_${ESSEX_KERNEL_LIB}

include ../../make.inc/rules
-include ${ESSEX_BASE}/make.inc/Makefile.${ESSEX_PLAT}

install_${ESSEX_KERNEL_LIB}:
	+cd ${ESSEX_KERNEL_LIB}; make install
