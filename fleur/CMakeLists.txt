cmake_minimum_required(VERSION 2.8)
project(Fleur)
enable_language(Fortran)

#add_definitions(-DDISPLAY_INFO)
##### Custom compile-time definitions explained:
# DISPLAY_INFO        Display lots of info during the runtime on stdout
# DUMP_DATA           REQUIRED for all data dumping, if any of the other flags are set,
#                       they will most likely not work without it (i.e. compile-time error)
# DUMP_DATA_OLD       enable data dumping using .bin files, works on simple systems
#                       with only a single atom type
# DUMP_DATA_HDF       new data dumping system to a file called 'fleur_dump.h5' using
#                       the HDF5 library
# DUMP_DATA_MT_HSSPHN_ONLY     Dump only the contribution from hssphn.F to
#                               Hamiltonian and Overlap matrix
# DUMP_ENABLE_SPHERICAL        Enable option to skip nonspherical part: If a file called
#                              .SPHERICAL_ONLY is found in an experiment directory, the
#                              matrices generated in hssphn.F only contain the spherical parts
# DUMP_SINGLE_KPT              Stop computation after the 1st k-point, ignoring if the kpts-file
#                              contains others
#
## Useful grep to find all DUMP_DATA.* uses:
# $ grep -n -P '(?<=#ifdef)\s+\bDUMP_DATA' -r . -A1
## All different dump flags in the code:
# $ grep -h -P '(?<=#ifdef)\s+\bDUMP_.*' -r . -o | sort | uniq
option(DUMP_RELEASE "Create data dumps when compiling in release mode" OFF)
option(DUMP_BIN "Dump internal data to .bin files for use in the C++ code" OFF)
option(DUMP_HDF "Create a single HDF5 file with all simulation input for the C++ code" ON)
option(DUMP_SPHERICAL "Enable skipping of nonspherical part (see CMakeLists.txt)" OFF)
option(DUMP_MT_ONLY "Dump the Muffin-Tin part of Hamiltonian and Overlap" ON)

# Fleur definitions from the original Makefile
add_definitions(-DCPP_DOUBLE
                -DCPP_APC
                -DCPP_F90)

# Adds data dump definitions described above if called
function(configure_data_dump)
    add_definitions(-DDUMP_DATA)
    if(DUMP_BIN)
        add_definitions(-DDUMP_DATA_OLD)
    elseif(DUMP_HDF)
        add_definitions(-DDUMP_DATA_HDF)
    endif()
    if(DUMP_SPHERICAL)
        add_definitions(-DDUMP_ENABLE_SPHERICAL)
    endif()
    if(DUMP_MT_ONLY)
        add_definitions(-DDUMP_DATA_MT_HSSPHN_ONLY)
    endif()
endfunction(configure_data_dump)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    configure_data_dump()
elseif(CMAKE_BUILD_TYPE MATCHES Release AND DUMP_RELEASE)
    configure_data_dump()
endif()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
#find_library(INTEL_IFPORT
#ifport
#PATHS ${INTEL_lib_root})
find_library(INTEL_IFCOREMT
    svml
    PATHS ${INTEL_lib_root})
find_library(INTEL_INTLC
    intlc
    PATHS ${INTEL_lib_root})

# TODO: find out what exactly those other libs are
# [JURECA] Removed IFPORT
#${INTEL_OMP} ${INTEL_IMF} ${INTEL_IFPORT} 
set(INTEL_MISC_LIBS 
    ${INTEL_OMP} ${INTEL_IMF}
    ${INTEL_IFCOREMT} ${INTEL_INTLC})

# HDF5 support, which needs MPI
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
#set(HDF5_INCLUDE_DIRS "/homeb/slai/slai15/hybridtcon/libs/hdf5-1.8.17/hdf5/include/")
#set(HDF5_Fortran_LIBRARIES "/homeb/slai/slai15/hybridtcon/libs/hdf5-1.8.17/hdf5/lib/libhdf5_fortran.so")
#set(HDF5_Fortran_HL_LIBRARIES "/homeb/slai/slai15/hybridtcon/libs/hdf5-1.8.17/hdf5/lib/libhdf5hl_fortran.so")
#set(HDF5_DEFINITIONS "/homeb/slai/slai15/hybridtcon/libs/hdf5-1.8.17/hdf5/include/H5public.h")
#set(HDF5_LIBRARIES ${HDF5_Fortran_LIBRARIES})
include_directories(${MPI_INCLUDE_PATH}
                    ${HDF5_INCLUDE_DIRS})
add_definitions(${HDF5_DEFINITIONS})
                
set(fleur_LIBS
        ${MKL_LIBS}
        ${MPI_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${HDF5_Fortran_HL_LIBRARIES})

    #message(WARNING ${HDF5_LIBRARIES})
    #message(WARNING ${HDF5_Fortran_HL_LIBRARIES})
    #message(WARNING ${HDF5_DEFINITIONS})


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)



# ifort  
if(Fortran_COMPILER_NAME MATCHES "ifort.*")
    message( "Chose ifort compiler, compiler name = ${Fortran_COMPILER_NAME}")
    # [JURECA] Added lstdc++
    set(CMAKE_Fortran_FLAGS 
    "-g -fopenmp -mkl -warn nousage -assume byterecl -r8 -W0 -fpp -heap-arrays -lstdc++")
    set(fleur_LIBS
            ${fleur_LIBS}
            ${INTEL_MISC_LIBS})
        #"-g -fopenmp -std90 -mkl -warn nousage -assume byterecl -r8 -W0 -fpp")
elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message( "Chose gfortran compiler, compiler name = ${Fortran_COMPILER_NAME}")
    set(CMAKE_Fortran_FLAGS
    "-g -fopenmp -fdefault-real-8 -fdefault-double-8 -m64 -cpp")
endif(Fortran_COMPILER_NAME MATCHES "ifort.*")

set(CMAKE_Fortran_FLAGS_RELEASE " -O4")
set(CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g")


set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_subdirectory(src)


#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur 
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)
