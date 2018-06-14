# Determine load balancing between GPUs and CPUs
About...

# Current state
...

# Building 

Requirements: 

 * Intel compiler icc with MKL
 * CUDA
 * C++11 compliant C++ compiler. GCC 4.8.2 will work
 * make.
 
On the RWTH cluster, use

    module load cuda

to get the required environment.

Buuilding is usually done as: 

 cd <root to ./testing/test_ratio>
 make mkl=1

The command generates the ./obj folder and builds <compute_ratio> programme.

# Running the code

