# High-performance implementation of FLAPW

This project aims to produce an implementation of 
the [FLAPW][] method in [DFT][]. Goals are:

* high performance: use (and research) how the mathematical equations
    can be implemented in a way that takes advantage of the architecture of 
    state-of-the-art multi-core and supercomputers.

* modularity and extensibility: follow modern programming paradigms,
    provide an easy-to-understand and tested code base

* portability: the code should be able to run on a laptop as well as a big cluster
    machine. 


[FLAPW]: http://webarchiv.fz-juelich.de/nic-series//volume31/bluegel.pdf 
    "Full-potential linearized augmented plane wave"
[DFT]: http://en.wikipedia.org/wiki/Density_functional_theory
    "Density functional theory"


# Current state

Generation of the H and S matrices is implemented using BLAS calls.  It uses
precomputed input froma specifically modified version of Fleur. The input data
includes the T-matrices, data about the system (atoms, types, cutoffs, energy
parameters), the k-point mesh the radial mesh and also the potential. 
For an exhaustive list, execute
    h5dump --onlyattr fleur_dump.h5 
on an input file produced by the modified Fleur version, as HDF5 is self-documenting.

The code takes these input quantities to generate the A and B tensors which are
then used in a high-performance BLAS kernel to generate H and S.
The relative error with regard to Fleur is on the order 10^-14, i.e. near machine accuracy.

# Building 

Requirements: 

 * Intel compiler ifort with MKL (for the wrapped Fortran code)
 * HDF5
 * C++11 compliant C++ compiler. GCC 4.8.2 will work, as will clang.
 
On the RWTH cluster, use

    module unload gcc
    module load gcc/4.8
    module load cmake
    module load LIBRARIES hdf5

to get the required environment.

CMake builds are usually done in an out-of-tree way: 

 cd <root directory of this repo>
 mkdir build && cd build
 # the next step 
 INTELROOT=/opt/intel/Compiler/default/compiler/ FC=ifort CXX=g++ CXXFLAGS=-O2 cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
 make -j16
 make install

This will configure, compile and install the code: In release mode, with prefix directory .local 
in the home folder.

## How to cite

If you use **HybridHSDLA** in your research, please cite the paper:

* Davor Davidović, Diego Fabregat-Traver, Markus Höhnerbach, Edoardo di Napoli (2018). *Accelerating the computation of FLAPW methods on heterogeneous architectures*, [[arXiv:1712.07206](https://arxiv.org/abs/1712.07206)].