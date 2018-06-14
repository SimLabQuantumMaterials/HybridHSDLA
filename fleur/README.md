What is this
============

A modified version of Fleur_v26 from http://flapw.de with extensions to dump the T matrices, potential
and other simulation parameters into a HDF5 file that can be read by the C++ implementation of the 
tensor contractions.
It also has been modified to use CMake as a build system so that builds can be done incremental and in parallel.


How to build
============


Requirements: 

 * Intel compiler ifort with MKL
 * HDF5, compiled with Fortran extensions for the Intel compiler
    -> correct .mod files, if HDF5 has been compiled with gfortran, it will not work with ifort

For the C++ part flapw_cpp:

 * C++11 compliant C++ compiler. GCC 4.8.2 will work, as will clang.


Use something like

INTELROOT=/opt/intel/Compiler/default/compiler/ HDF5_ROOT=~/.local FC=ifort CXX=clang++ CXXFLAGS=-O3 cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local ..

from a build directory. 

INTELROOT - should point at the installation location of the Intel compiler, i.e. have subfolders
            lib/ and mkl/

HDF5_ROOT - needed if the default installation that cmake can find does not contain the proper Fortran extensions.
            Specify the install prefix of the HDF5 installation here.
