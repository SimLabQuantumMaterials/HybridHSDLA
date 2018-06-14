#!/usr/bin/env bash

# The following environment variables need to be
# properly setup for this script to work:
# - INTELROOT
# - HDF5_ROOT

LIBS=NO
FLEUR=NO
CPU=NO
HYBRID=NO
HS_CUDA=NO
HS_MIC=NO

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        LIBS)
        LIBS=YES
        ;;
        FLEUR)
        FLEUR=YES
        ;;
        CPU)
        CPU=YES
        ;;
        HYBRID)
        HYBRID=YES
        ;;
        HS_CUDA)
        HS_CUDA=YES
        ;;
        HS_MIC)
        HS_MIC=YES
        ;;
        *)
        echo "[ERROR] Wrong option $key"
        echo "Usage $0 LIBS|FLEUR|CPU|HYBRID|HS_CUDA|HS_MIC"
        exit 1
        ;;
    esac
    shift # to the next argument
done

if [ $LIBS = YES ]
then
    ###############
    # HDF5
    ###############
    [ -e libs ] || mkdir libs/
 
    procs=`grep -c ^processor /proc/cpuinfo`
    hdf5_version="1.10.1"

    cd libs
    # Download hdf5 1.10 (or newer stable version)
    hdf5_base=`echo ${hdf5_version} | sed -r 's/.{2}$//'`
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_base}/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.bz2
    # Unpack
    tar -vjxf hdf5-${hdf5_version}.tar.bz2
    # Compile it
    cd hdf5*
    #FC=ifort ./configure --enable-fortran --enable-fortran2003
    make -j $procs
    make install
    cd ../..

    ###############
    # Armadillo
    ###############
    cd flapw_opt/external
    arma_version="8.500.1"
    wget http://sourceforge.net/projects/arma/files/armadillo-${arma_version}.tar.gz
    tar -xzf armadillo-${arma_version}.tar.gz
    mv armadillo-${arma_version} armadillo
    cd ../..
fi    


########
# FLEUR
########
if [ $FLEUR = YES ]
then
    cd fleur
    [ -e build ] || mkdir build 
    cd build

    CC=icc FC=ifort \
    CXX=icpc CXXFLAGS=-O3 \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local .. -DDUMP_RELEASE=ON

    make -j 8

    cd ../../
fi

############
# CPU stripped version
############
if [ $CPU = YES ]
then
    # GCC >= 4.8 !!
    cd flapw_opt
    [ -e build_cpu ] || mkdir build_cpu 
    cd build_cpu

    CC=icc FC=ifort CXX=icpc CXXFLAGS=-O3 \
    cmake -DTARGET_IMPLEMENTATION=CPU_STRIPPED -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local VERBOSE=1 ..

    make -j 8

    cd ../..
fi

##########
# HYBRID
##########
if [ $HYBRID = YES ]
then
    cd flapw_opt
    [ -e build_hybrid ] || mkdir build_hybrid 
    cd build_hybrid

    CC=icc CFLAGS=-O3 FC=ifort \
    CXX=icpc CXXFLAGS=-O3 \
    cmake -DTARGET_IMPLEMENTATION=HYBRID -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local VERBOSE=1 ..

    make -j 8

    cd ../..
fi

############
# hs_cuda
############
if [ $HS_CUDA = YES ]
then
    cd flapw_opt
    [ -e build_hs_cuda ] || mkdir build_hs_cuda
    cd build_hs_cuda

    CC=icc FC=ifort \
    CXX=icpc CXXFLAGS=-O3 \
    cmake -DTARGET_IMPLEMENTATION=HS_CUDA -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local VERBOSE=1 ..

    make -j 8

    cd ../..
fi

############
# hs_mic
############
if [ $HS_MIC = YES ]
then
    cd flapw_opt
    [ -e build_hs_mic ] || mkdir build_hs_mic
    cd build_hs_mic

    CC=icc FC=ifort \
    CXX=icpc CXXFLAGS=-O3 \
    cmake -DTARGET_IMPLEMENTATION=HS_MIC -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local VERBOSE=1 ..

    make -j 8

    cd ../..
fi
