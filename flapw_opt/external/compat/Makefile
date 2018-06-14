DEFAULT_FC = ifort
DEFAULT_FC_FILE = Make.$(DEFAULT_FC).inc

ifeq ($(FC), ifort)
	include Make.ifort.inc
else ifeq ($(FC), gfortran)
	include Make.gfortran.inc
else
	FC = $(DEFAULT_FC)
	include $(DEFAULT_FC_FILE)
endif

CPPDEFS = -DCPP_DOUBLE -DCPP_APC

# source files from original fleur
VPATH = fleur

FCFLAGS += $(CPPDEFS)

.PHONY: all copy clean

all: rad_test libflerp.so

copy: libflerp.so
	cp $^ ../

clean:
	rm -f *.o *.mod *.so rad_test libflerp.so

# a bit clunky, but who's got time to write pretty rules?
# recompiles all objects as position-independent code (PIC)
# so that they can be linked into the shared library

libflerp.so: libflerp.pic.o radfun.pic.o \
	constants.pic.o radsra.pic.o intgr.pic.o radsrd.pic.o \
	ylm4.pic.o dsphbs.pic.o sphbes.pic.o

# ugly: .mod files for objects shared with rad_test exist, since it makes no difference
# for .mod if the object is PIC or not. so, explicitly add those from libflerp.so to libflerp.pic.so

libflerp.pic.o: ylm4.pic.o dsphbs.pic.o sphbes.pic.o

rad_test: rad_test.o radfun.o constants.o radsra.o intgr.o radsrd.o binmat.o

radfun.o: constants.o radsra.o intgr.o radsrd.o




rad_test.o: radfun.o binmat.o

%.so: 
	$(FC) $(FCFLAGS) -shared -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.pic.o: %.f90
		$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

%.pic.o: %.f
		$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

%.pic.o: %.F
		$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F
	$(FC) $(FCFLAGS) -o $@ -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)


radsra.o: intgr.o

