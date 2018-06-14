# Directories
CUDA_HOME ?= /rwthfs/rz/SW/cuda/7.5.18/RHEL_7/cuda
MAGMA_HOME ?= /usr/local/magma-1.6.1
LAPACK ?= /usr/local/lapack-3.5.0
BLAS ?= /usr/local/OpenBLAS/0.2.18/multiThread/
CUDA_LIB ?= $(CUDA_HOME)/lib64
MAGMA_LIB ?= $(MAGMA_HOME)/lib
CUDA_INCLUDE ?= $(CUDA_HOME)/include
MAGMA_INCLUDE ?= $(MAGMA_HOME)/include
BLAS_INCLUDE ?= $(BLAS)/include
BLAS_LIB ?= $(BLAS)/lib

F77 	= ifort
CC      = icc
CLINKER = icpc
CFLAGS  = -Wall -fopenmp
DBG 	= 

ifeq ($(dbg),1)
	DBG += -g -O0 -fopenmp
else
	FFLAGS += -O3
	CFLAGS += -O3
endif

ifeq ($(timer),1)
	CFLAGS += -DTIMER
endif

ifeq ($(mkl),1)
    CFLAGS += -DMKL
endif

LIBDIRS = -L/usr/lib64 -L$(LAPACK) -L$(BLAS_LIB) -L$(CUDA_LIB)
INCDIRS = -I$(CUDA_INCLUDE) -I$(BLAS_INCLUDE)
ODIR = ./obj
SDIR = .
RDIR = .

CUDA_LIBS = -lcudart -lcublas
OPENBLAS_LIBS = -lgfortran -lopenblas -llapack
OTHER_LIBS = -lm -lpthread
LIBLIST = -lcudart -lcublas -lgfortran -lopenblas -llapack -lm -lpthread
MKL_LIBS = -mkl

_EXECS = compute_ratio
EXECS = $(patsubst %,$(RDIR)/%,$(_EXECS))

_OBJ_GPU = compute_ratio.o kernels.o
_OBJ_OTH = timing.o

OBJ_GPU = $(patsubst %,$(ODIR)/%,$(_OBJ_GPU))
OBJ_OTH = $(patsubst %,$(ODIR)/%,$(_OBJ_OTH))

dir_create=@mkdir -p $(@D)

default: $(EXECS)
all: $(EXECS) 

$(RDIR)/compute_ratio: $(ODIR)/compute_ratio.o $(OBJ_GPU) $(OBJ_OTH)
	$(dir_create)
	$(CC) $(CFLAGS) $(LIBDIRS) -o $@ $^ $(CUDA_LIBS) $(MKL_LIBS) $(OTHER_LIBS)

$(ODIR)/%.o: $(SDIR)/%.c
	echo Compiling C files...
	$(dir_create)
	$(CC) $(CFLAGS) $(DBG) $(INCDIRS) -c -o $@ $<	

$(ODIR)/%.o: $(SDIR)/%.cpp
	echo Compiling C++ files...
	$(dir_create)
	$(CC) $(CFLAGS) $(DBG) $(INCDIRS) -c -o $@ $<	

.PHONY: clean

clean:
	rm -f $(EXECS) $(ODIR)/*.o 

clean_all:
	rm -f $(EXECS) *.o
