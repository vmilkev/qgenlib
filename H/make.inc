#================================================================
#  hmat make include file.
#  version 0 releas 0.0
#  May 2018
#================================================================

SHELL = /bin/sh
SIZE=

#===============================================================

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
OS = linux
else
OS = _Win64
endif


CXX= icc
MATRLIB = ${MY_PROJ}UTILITIES/QGMATHLIB/qgenlib/matrixclass/mathlib/lib/matrix.a


# debug configuration:
#CFLAGS = -g -O0 -qopenmp -debug parallel -fp-stack-check -ftrapuv -fargument-alias -check=uninit -check=stack -check=conversions -traceback -Wall -D $(OS) -MMD -std=c++17


# release configuration:
# NOTE: options assumes 64 based architecture, 32 bit integer interface layer

CppFLAGS = -O3 -qopenmp -mcmodel=large -I${MATRLIB} -I${MKLROOT}/include -std=c++17 -traceback -Wall -qoverride-limits -D $(OS) -MMD

CppLINK = -cxxlib -static-libgcc -static-libstdc++ -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
# CppLINK = -cxxlib -qopenmp -L${MKLROOT}/lib/intel64_lin -lmkl_def -lmkl_avx2 -lm -ldl

# Static linking:
# LFLAGS = ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

#================================================================
#
#  Set directory names
SRC = src/

VPATH = $(SRC)
vpath %.o  lib/

BIN= bin
MYLIB= lib/
F95ROOT=/opt/intel/current
