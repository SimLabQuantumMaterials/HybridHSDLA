/*
 * zherkx_hyb.c
 *
 *  Created on: Mar 22, 2017
 *      Author: Davor Davidovic
 *      Rudjer Boskovic Institute
 *
 *      The algorithms performs one of the hermitian rank k operations
 *
 *      C := alpha*A*A**H + beta*C,
 *
 *      or
 *
 *      C := alpha*A**H*A + beta*C,
 *
 *      where  alpha and beta are real scalars, C is an n x n hermitian
 *      matrix and  A is an n x k matrix in the first case and a k x n
 *      matrix in the second case.
 *    
 *      The work is divided by the GPUs and CPUs as follows:
 *      The leading principal submatrix of C is executed on the GPUs (cublasXtZherk) while the rest of
 *      the matrix C is updated on the CPUs (multithreaded/openmp).
 *      
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cublasXt.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <omp.h>
#include "utils.h"

#ifdef MKL
  #include <mkl_service.h>
  #include <mkl.h>
#else
  #include <cblas.h>
#endif

#define max(a,b) (((a)>(b))?(a):(b))


int zherk_hyb( cublasXtHandle_t handle, 
               char uplo, char trans, 
               int n, int k, 
               double alpha,
               double complex *A, int lda,  
               double beta,
               double complex *C, int ldc, 
               int nb, double ratio, int nthreads )
{
/* 
  Arguments
  ==========

  See: http://docs.nvidia.com/cuda/cublas/#cublas-lt-t-gt-herkx

RATIO    - On entry, RATIO specifies the ratio of matrix C that
           will be executed on the GPUs. If the RATIO = 0,
	         work will be equally split between GPUs and CPUs.

NTHREADS - Number of threads that will be used in hybrid exection.
	         If NTHREADS == 0 then the total number of threads is set
	         to the max number of threads returned by MKL/OpenBLAS
*/


/* .. CUBLAS types .. */
  cublasOperation_t cu_trans;
  cublasOperation_t cu_transa;
  cublasFillMode_t  cu_uplo;
  cublasStatus_t    status;

/* .. CUBLASXT variables */
//  int *devices = NULL;
//  cublasXtHandle_t handle;

/* .. Local Scalars .. */
  int  info;
  char *mesg;
  char transa, transb;
  int  GPU_n, CPU_n;
  int  n_cpus, n_gpus;
  int  my_id;
  int  offset;

/* .. Parameters .. */
  static double complex ONEZ = 1.0 + 0.0*I, ZEROZ = 0.0 + 0.0*I;


/* =======  Executable part ======= */


/*	Test the input parameters.	*/
  
  info = 0;
  if ( n < 0 )
    info = 4;
  else if ( k < 0 )
    info = 5;
  else if ( trans == 'N' && lda < max(1,n) )
    info = 8;
  else if ( trans != 'N' && lda < max(1,k) )
    info = 8;
  else if ( ldc < max(1,n) )
    info = 11;
  else if ( nb < 0 )
    info = 12;
  else if ( ratio < 0.0 )
    info = 13;
  else if ( nthreads < 0 )
    info = 14;
  if ( info != 0 ){
    mesg = "ZHERK_HYB";
    xerbla(mesg, &info, 9);
    return -1*info;
  }

/*	Quick return if possible	*/

  if ( ( k == 0 ) || ( ( alpha == ZEROZ ) && ( beta == ONEZ ) ) )
    return info;
  

/*  If alpha == zero scale C with beta on CPU  */

  if( alpha == ZEROZ ){
    if( trans == 'N'){
      transa = 'N';
    }else{
      transa = 'C';
    }
    
    zherk(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);

  return info;
  }
  
/*  .. Preparation part  .. */

/*  Get the number of GPU devices  */
  LOG_CUDA_STATUS( cudaGetDeviceCount(&n_gpus), "cudaGetDeviceCount" );

/*  Get number of threads, if number of threads not given, get max number of threads */
  if( nthreads == 0){
    n_cpus = omp_get_max_threads();
  }else{
    n_cpus = nthreads;
  }

  printf("ZHERK_HYB is running on %d CPUs and %d GPUs\n", n_cpus, n_gpus);

/*  If ratio is not set, set it to 0.5 (equally distribute between GPU and the CPU  */
  if( ratio == 0.0 )
    ratio = 0.5;

/*  Compute the number of rows of the matrix C whose update will be offloaded to the GPU/CPU  */
  GPU_n = n/ratio;
  CPU_n = n - GPU_n;

  printf("ZHERK_HYB: GPU part (%d x %d)\n", GPU_n, GPU_n);
  printf("ZHERK_HYB: CPU part (%d x %d)\n", CPU_n, n);

/*  Set openmp and mkl variables for nested parallelism  */
#ifdef MKL
  mkl_set_dynamic(0);
#endif
  omp_set_num_threads(2);
  omp_set_nested(1);


/*  Start parallel section. Two threads, one for GPU, one for CPU  */

#pragma omp parallel default(shared) private(my_id, cu_uplo, cu_trans)
{

  my_id = omp_get_thread_num();

/*  Thread id == 0 manages GPU execution  */
  if( my_id == 0){

/*  Convert to CUBLAS types for UPLO and TRANS */
    cu_uplo  = uplo  == 'L' ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
    cu_trans = trans == 'N' ? CUBLAS_OP_N : 
               trans == 'T' ? CUBLAS_OP_T : CUBLAS_OP_C;

/*  Call CUBLASXT routine  */
    //status = cublasXtZherk(handle, cu_uplo, cu_trans, GPU_n, k, &alpha, (cuDoubleComplex *)A, lda, &beta, (cuDoubleComplex *)C, ldc);
    LOG_CUBLAS_STATUS( cublasXtZherk(handle, cu_uplo, cu_trans, GPU_n, k, &alpha, (cuDoubleComplex *)A, lda, &beta, (cuDoubleComplex *)C, ldc), "cublasXtZherk" );

  }else{
/*  I'm thread 1 and manage execution on CPUs  */   

/*  First compute the diagonal block (call to zherk) */
    if( trans == 'N' ){
      offset = GPU_n;
      transa = 'N';
    }else{
      offset = GPU_n * lda;
      transa = 'C';
    }

/*  Set number of threads for multithread blas routins (-1 thread that is reserverd for the GPU management) */
#ifdef MKL
    mkl_set_num_threads(n_cpus-1);
#else
    openblas_set_num_threads(n_cpus-1);
#endif

    zherk(&uplo, &transa, &CPU_n, &k, &alpha, &A[offset], &lda, &beta, &C[GPU_n + GPU_n*ldc], &ldc);

/*  Then compute the last row-panel  */
    if( trans == 'N' ){
      offset = GPU_n;
      transa = 'N';
      transb = 'C';
    }else{
      offset = GPU_n * lda;
      transa = 'C';
      transb = 'N';
    }

    double complex alphaComplex, betaComplex;

    alphaComplex = alpha;
    betaComplex = beta;
    zgemm(&transa, &transb, &CPU_n, &GPU_n, &k, &alphaComplex, &A[offset], &lda, A, &lda, &betaComplex, &C[GPU_n], &ldc);

  }
}

  /* Return the number of threads for blas operations that was before this function call */
#ifdef MKL
  mkl_set_num_threads(n_cpus);
#else
  openblas_set_num_threads(n_cpus);
#endif
  omp_set_num_threads(n_cpus);
  
  return info;
}
