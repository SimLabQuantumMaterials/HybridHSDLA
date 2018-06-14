Heteregeneous implementation of ZGEMM, ZHERK and ZHER2K.

Consists of the files
- source.h
- source.cpp
- sink.c

Supports
- NVIDIA GPUS using CUDA (tested w/ version 7.5.18) 
- Intel Xeon Phi accelerator cards using hStreams (tested w/ commit 06dd00351ac122290980b760d0b3cd38ca3c4e62, path to installation needs to be set using $HETERO)
- Some BLAS on the host (Makefile uses MKL by default)

Can be tuned using environment variables
- to the max. block size: HSBS=bs
- to the n-fold buffering: HSNB=n
- to if the host participates: HSHC=0/1
- to the depth of initial warm-up using smaller blocks: HSSPLIT=...
- to the number of devices to use: HSND=...

Possible future tasks
- Generalize to other BLAS operations (very straightforward, but not necessary right now).
- Benchmark individual operations.
- Limit advance() to one unpack per call.
- Eliminate packing in GPU case.
- Add thresholds to switch into host-BLAS if dims too small
- Tune to find better defaults for params.
- Add functions to register and unregister matrices, i.e. they will only be written to using library functions, so we need to wait for calcs. to complete not n times, but just once per assembly.
- Limit to one unpack operation per issue call.
- Improve block size selection.
- (maybe) have fun with LD_PRELOAD and stuff s.th. the library works for programs it is not compiled for.
