#FC = gfortran
FCFLAGS = -fdefault-real-8 -g
LDFLAGS= -L$(MKLROOT)/lib/intel64
#LDLIBS = -lmkl_gf_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
LDLIBS = -lmkl_rt -lpthread -lm
