/*
 * zherkx_hyb.c
 *
 *  Created on: July 20, 2017
 *      Author: Davor Davidovic
 *      Rudjer Boskovic Institute
 *
 *      ZHER2K routine performs one of the Hermitian rank 2k operations
 *
 *      C := alpha*A*B**H + conjg(alpha)*B*A**H + beta*C,
 *
 *      or
 *
 *      C := alpha*A**H*B + conjg(alpha)*B**H*A + beta*C,         
 *
 *      where alpha and beta are scalars with beta real, C is an n-by-n
 *      Hermitian matrix stored in lower or upper mode, and A and B are 
 *      n-by-k matrice in the first case and k-by-n matrice in the second
 *      case.
 *    
 *      The work is divided between GPUs and CPUs as follows:
 *      The leading principal submatrix of C is executed on the GPUs (cublasXtZher2k) 
 *      while the rest of the matrix is updated on the CPUs (multithreaded mkl/openblas).
 *       ______________   __________________   _____________     __________________   _____________
 *      |      |      |   |                |   |     |     |     |                |   |     |     |    
 *      | C_00 | C_01 |   |     A_0        |   |     |     |     |       B_0      |   |     |     |
 *      |______|______|   |________________|   |     |     |     |________________|   |     |     |
 *      |      |      | = |                | * |     |     |     |                | * |     |     |
 *      | C_10 | C_11 |   |     A_1        |   | B_0 | B_1 |  +  |       B_1      |   | A_0 | A_1 |
 *      |______|______|   |________________|   |     |     |     |________________|   |     |     |
 *                                             |     |     |                          |     |     |
 *                                             |     |     |                          |     |     |
 *                                             |_____|_____|                          |_____|_____|
 *
 *             C         =        A          *       B**H     +          B          *       A**H
 *
 *      In case, C is stored in lower mode and A and B are n-by-k case, then C is computed as follows
 *
 *      C_00 = A_0 * B_0**H  +  B_0 * A_0**H  --> 1 call to zher2k (GPU) 
 *      C_11 = A_1 * B_1**H  +  B_1 * A_1**H  --> 1 call to zherk2k (CPU)
 *      C_10 = A_1 * B_0**H  +  B_1 * A_0**H  --> 2 calls to zgemm (CPU)
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

#define TRUE 1
#define FALSE 0
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

int zher2k_hyb( cublasXtHandle_t handle, 
               char uplo, char trans, 
               int n, int k, 
               double complex alpha,
               double complex *A, int lda,
               double complex *B, int ldb,  
               double beta,
               double complex *C, int ldc, 
               int nb, double ratio, int nthreads )
{
/* 
  Arguments
  ==========

  See: http://docs.nvidia.com/cuda/cublas/#cublas-lt-t-gt-herkx

RATIO    - On entry, RATIO specifies the ratio of matrix C that
           will be executed on the GPUs. If the RATIO = 0.5,
	   the work will be equally split between the GPUs and CPUs.

NTHREADS - Number of threads that will be used in hybrid exection.
	   If NTHREADS == 0 then the total number of threads is set
	   to the max number of threads returned by MKL/OpenBLAS
*/


/* .. CUBLAS types .. */
  cublasOperation_t cu_trans;
  cublasFillMode_t  cu_uplo;
  cublasStatus_t    status;

/* .. Local Scalars .. */
  int  info;
  char *mesg;
  char transa, transb;
  int  GPU_n, CPU_n;
  int  n_cpus, n_gpus;
  int  my_id;
  int  offset_a, offset_b, offset_c;

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
  else if ( trans == 'N' && ldb < max(1,n) )
    info = 10;
  else if ( trans != 'N' && lda < max(1,k) )
    info = 10;
  else if ( ldc < max(1,n) )
    info = 13;
  else if ( nb < 0 )
    info = 14;
  else if ( ratio < 0.0 )
    info = 15;
  else if ( nthreads < 0 )
    info = 16;
  if ( info != 0 ){
    mesg = "ZHER2K_HYB";
    xerbla(mesg, &info, 10);
    return -1*info;
  }

/*	Quick return if possible	*/

  if ( ( n == 0 ) || ( ( ( alpha == ZEROZ ) || ( k == 0 ) ) && ( beta == 1.0) ) )
    return info;
  

/*  If alpha == zero scale C with beta on CPU  */

  if( alpha == ZEROZ ){
    
    zher2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);

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

  printf("ZHER2K_HYB is running on %d CPUs and %d GPUs\n", n_cpus, n_gpus);

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


/*  Start parallel section. Two master threads, one for GPU, one for CPU  */

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
//    status = cublasXtZher2k(handle, cu_uplo, cu_trans, GPU_n, k, (const cuDoubleComplex*)&alpha, (cuDoubleComplex *)A, lda, (cuDoubleComplex *)B, ldb, &beta, (cuDoubleComplex *)C, ldc);
    LOG_CUBLAS_STATUS( cublasXtZher2k(handle, cu_uplo, cu_trans, GPU_n, k, (const cuDoubleComplex*)&alpha, (cuDoubleComplex *)A, lda, (cuDoubleComplex *)B, ldb, &beta, (cuDoubleComplex *)C, ldc), "cublasXtZherk2k" );

  }else{
/*  I'm thread 1 and manage execution on CPUs  */   

/*  First compute the diagonal block (call to zher2k) */
    if( trans == 'N' ){
      offset_a = offset_b = GPU_n;
      transa = 'N';
      transb = 'C';
    }else{
      offset_a = GPU_n * lda;
      offset_b = GPU_n * ldb;
      transa = 'C';
      transb = 'N';
    }
    offset_c = GPU_n + GPU_n * ldc;

/*  Set number of threads for multithread blas routins (-1 thread that is reserverd for the GPU management) */
#ifdef MKL
    mkl_set_num_threads(n_cpus-1);
#else
    openblas_set_num_threads(n_cpus-1);
#endif

    zher2k(&uplo, &transa, &CPU_n, &k, &alpha, &A[offset_a], &lda, &B[offset_b], &ldb, &beta, &C[offset_c], &ldc);

/*  Then compute the last row-panel (C_10) - 2 zgemm operations  
 *
 *  C_10 = A_1 * B_0**H + B_1 * A_0**H + beta * C_10 (trans = 'N')
 *  C_10 = A_1**H * B_0 + B_1**H * A_0 + beta * C_10 (trans = 'C')
 *
 *  The computation is divided in two consecutive zgemm calls
 *  C_10 = B_1 * A_0**H + beta * C_10 (trans = 'N') and C_10 = B_1**H * A_0 + beta * C_10 (trans = 'C') (1st zgemm)
 *  C_10 = A_1 * B_0**H + 1 + C_10 (trans = 'N') and C_10 = A_1**H * B_0 + 1 * C_10 (trans = 'C') (2nd zgemm)
 *
 */
  
    double complex alphaConj, betaComplex;
    alphaConj = conj(alpha);
    betaComplex = beta;

    if( uplo == 'L' || uplo == 'l' ){
	offset_c = GPU_n;

	/* First zgemm call */
    	zgemm(&transa, &transb, &CPU_n, &GPU_n, &k, &alphaConj, &B[offset_b], &ldb, A, &lda, &betaComplex, &C[offset_c], &ldc);
    	/* Second zgemm call */
    	zgemm(&transa, &transb, &CPU_n, &GPU_n, &k, &alpha, &A[offset_a], &lda, B, &ldb, &ONEZ, &C[offset_c], &ldc);

    }else{
	offset_c = GPU_n * ldc;

    	/* First zgemm call */
    	zgemm(&transa, &transb, &GPU_n, &CPU_n, &k, &alphaConj, B, &ldb, &A[offset_a], &lda, &betaComplex, &C[offset_c], &ldc);

	/* Second zgemm call */
    	zgemm(&transa, &transb, &GPU_n, &CPU_n, &k, &alpha, A, &lda, &B[offset_b], &ldb, &ONEZ, &C[offset_c], &ldc);
    }
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
