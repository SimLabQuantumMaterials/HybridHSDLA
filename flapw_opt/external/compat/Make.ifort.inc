#ifort -o build/debug/radfun.o -c -fopenmp -O4 -g -mkl -warn nousage -assume byterecl -r8 -W0 -fpp -O1 -g -DCPP_DOUBLE -DCPP_APC -DDUMP_DATA -Ibuild/debug -module build/debug build/debug/radfun.f
#FC = ifort
FCFLAGS = -fopenmp -O4 -g -mkl -warn nousage -assume byterecl -r8 -W0 -fpp -O1 
LDFLAGS= -L$(MKLROOT)/lib/intel64
LDLIBS = -lmkl_rt -lmkl_def -lpthread -lm
