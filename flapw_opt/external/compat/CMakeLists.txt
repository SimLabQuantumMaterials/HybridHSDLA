cmake_minimum_required(VERSION 2.8)
project(FLERP Fortran) # FLEur wRaPper... derp!

include_directories(${PROJECT_SOURCE_DIR}/fleur)
add_definitions("-DCPP_DOUBLE -DCPP_APC")

# Intel libraries

set(INTEL_root $ENV{INTELROOT})
set(MKL_lib_root ${INTEL_root}/mkl/lib/intel64;$ENV{MKLROOT}/lib/intel64)
set(INTEL_lib_root ${INTEL_root}/lib/intel64)

# MKL
find_library(MKL_RT 
    mkl_rt 
    PATHS ${MKL_lib_root})
#find_library(MKL_DEF 
#mkl_def 
#PATHS ${MKL_lib_root})
#set(MKL_LIBS ${MKL_RT} ${MKL_DEF} pthread m)
set(MKL_LIBS ${MKL_RT} pthread m)
# others (OMP, ...)
find_library(INTEL_OMP
    iomp5
    PATHS ${INTEL_lib_root})
find_library(INTEL_IMF
    imf
    PATHS ${INTEL_lib_root})
find_library(INTEL_IFPORT
    ifport
    PATHS ${INTEL_lib_root})
find_library(INTEL_IFCOREMT
    svml
    PATHS ${INTEL_lib_root})
find_library(INTEL_INTLC
    intlc
    PATHS ${INTEL_lib_root})

# TODO: find out what exactly those other libs are
#set(INTEL_MISC_LIBS 
#${INTEL_OMP} ${INTEL_IMF} ${INTEL_IFPORT} 
#${INTEL_IFCOREMT} ${INTEL_INTLC})
set(INTEL_MISC_LIBS 
    ${INTEL_OMP} ${INTEL_IMF}
    ${INTEL_IFCOREMT} ${INTEL_INTLC})
                

set(flerp_LIBS ${MKL_LIBS} ${INTEL_MISC_LIBS})

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  message("Chose gfortran compiler, compiler name = ${Fortran_COMPILER_NAME}")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -g")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   " -O1")
elseif(Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort  
  message( "Chose ifort compiler, compiler name = ${Fortran_COMPILER_NAME}")
  set(CMAKE_Fortran_FLAGS 
  "-g -fopenmp -g -mkl -warn nousage -assume byterecl -r8 -W0 -fpp")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O4")
  set(CMAKE_Fortran_FLAGS_DEBUG   " -O1 -g")
  #message("MKLROOT = $ENV{MKLROOT}")
endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")


set(flerp_SRC 
  libflerp.f90 
  fleur/radfun.f 
  fleur/constants.f
  fleur/radsra.f 
  fleur/intgr.F 
  fleur/radsrd.F 
  fleur/ylm4.f 
  fleur/dsphbs.f 
  fleur/sphbes.f)

add_library(flerp SHARED ${flerp_SRC})
target_link_libraries(flerp ${flerp_LIBS})

install(TARGETS flerp
        LIBRARY DESTINATION lib
        CONFIGURATIONS Release)
