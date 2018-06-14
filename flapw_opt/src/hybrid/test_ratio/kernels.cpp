/*
 * ratio_zherk.c
 *
 *  Created on: July 25, 2017
 *      Author: Davor Davidovic
 *      Rudjer Boskovic Institute
 *
 *      The rutine compares the performance between the CPU and GPU devices and returns the ratio (GPU_perf/CPU_perf),
 *      i.e. what is the speedup of GPU compared to () 
 * 
 *      The comparison is done by performing several BLAS routines on both CPU and GPU and 
 *      comparing their performance
 */

#include <stdio.h>
#include <complex>
#include <cublas_v2.h>
#include <math.h>
#include <cublasXt.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#ifdef MKL
  #include <mkl_service.h>
#else
  #include <cblas.h>
#endif

extern "C"{
    #include "timings.h"
    void zherk_(char*, char*, int*, int*, double*, std::complex<double>*, int*, double*, std::complex<double>*, int*);
    void zher2k_(char*, char*, int*, int*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*);
#ifdef MKL
    void zgemmt_(char*, char*, char*, int*, int*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*);
#endif
    void zgemm_(char*, char*, int*, int*, int*, double*, std::complex<double>*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*);
}

#define MAX_DEV 16
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

void ratio_zherk(cublasXtHandle_t handle, int min_n, int max_n, int max_k, int step, std::complex<double> *A, int lda, std::complex<double> *C, int ldc, double *ratio)
{

    int    ndevices, devices[MAX_DEV], i;
    std::complex<double> alpha = 1.66 - 0.43i;
    std::complex<double> beta = -3.14 + 1.97i;
    double alpha_d, beta_d;
    char   transa, transb, uplo;
    int    n, k, GPU_n, CPU_n, offset;
    int    start, end;
    double gpu_gflops, cpu_gflops;
    double exec_time, best_time;

    /*  CUBLAS */
    cublasOperation_t cu_trans;
    cublasFillMode_t cu_uplo;
    cublasStatus_t status;
    
    /*  Timing routines  */
    struct timeval t_start, t_end;
    long int gpu_t, cpu_t;

    /*  Get number of cuda devices */
    cudaGetDeviceCount( &ndevices );

    /*  Set devices */
    for( i = 0; i < ndevices; i++)
        devices[i] = i;

    /*  Set parameters to zherk calls */
    cu_trans = CUBLAS_OP_N;
    cu_uplo = CUBLAS_FILL_MODE_LOWER;
    transa = 'N';
    transb = 'C';
    uplo = 'L';
    start = min_n;
    end = max_n;
    alpha_d = alpha.real();
    beta_d = beta.real();

/*  Run a series of test with the increasing number of rows of C offloaded to the GPU  */

    printf("\n");
    printf("dev    GPU_N   CPU_N   K    Time(sec) GFLOPS\n");
    printf("%s\n", std::string(45,'-').c_str());
    best_time = 1.0e3;
    for( i = start; i <= end; i+=step){

        GPU_n = i;
        CPU_n = max_n - GPU_n;

        /* Compute GPU part */
        read_clock(&t_start);

    	status = cublasXtZherk(handle, cu_uplo, cu_trans, GPU_n, max_k, &alpha_d, (const cuDoubleComplex *)A, lda, &beta_d, (cuDoubleComplex *)C, ldc);
    	if(status != CUBLAS_STATUS_SUCCESS){
    		printf("Error in %s:%d CublasXtZherk error! %d\n", __FILE__, __LINE__, status);
    	}
        
        cudaThreadSynchronize();
        read_clock(&t_end);
        gpu_t = elapsed_time(&t_start, &t_end);
        gpu_gflops = (4.0 * (double)GPU_n * (double)max_k * ((double)GPU_n + 1))/(gpu_t*1e3);

	/* Compute CPU part */
    	read_clock(&t_start);

	    offset = GPU_n;

        zherk_(&uplo, &transa, &CPU_n, &max_k, &alpha_d, &A[offset], &lda, &beta_d, &C[GPU_n + GPU_n*ldc], &ldc);
        zgemm_(&transa, &transb, &CPU_n, &GPU_n, &max_k, &alpha_d, &A[offset], &lda, A, &lda, &beta_d, &C[GPU_n], &ldc);

    	read_clock(&t_end);
    	cpu_t = elapsed_time(&t_start, &t_end);
        cpu_gflops = (4.0 * CPU_n * max_k * (CPU_n + 1 + 2*GPU_n))/(cpu_t*1e3);

    	printf("    %6d   %6d %6d\n", GPU_n, CPU_n, max_k);
        printf("GPU                        %7.4f   %7.4f\n", gpu_t*1e-6, gpu_gflops);
        printf("CPU                        %7.4f   %7.4f\n", cpu_t*1e-6, cpu_gflops);

        exec_time = max(gpu_t, cpu_t)*1e-6;
//	    error = fabs(gpu_t - cpu_t)/cpu_t;

        if( exec_time < best_time ){
	        *ratio = (double)max_n / (double)GPU_n;
	        best_time = exec_time;
	    }
    }
    
    return;
}


void ratio_zherkx( cublasXtHandle_t handle, int min_n, int max_n, int max_k, int step, std::complex<double> *A, int lda, std::complex<double> *B, int ldb, std::complex<double> *C, int ldc, double *ratio )
{
    int    ndevices, devices[MAX_DEV], i;
    std::complex<double> alpha = 1.66 - 0.43i;
    std::complex<double> beta = -3.14 + 1.97i;
    double alpha_d, beta_d;
    char   transa, transb, uplo;
    int    n, k, GPU_n, CPU_n, offset;
    int    start, end;
    double gpu_gflops, cpu_gflops;
    double exec_time, best_time;

    /*  CUBLAS */
    cublasOperation_t cu_trans;
    cublasFillMode_t cu_uplo;
    cublasStatus_t status;
    
    /*  Timing routines  */
    struct timeval t_start, t_end;
    long int gpu_t, cpu_t;

    /*  Get number of cuda devices */
    cudaGetDeviceCount( &ndevices );

    /*  Set devices */
    for( i = 0; i < ndevices; i++)
        devices[i] = i;

    /*  Set parameters to zherk calls */
    cu_trans = CUBLAS_OP_N;
    cu_uplo = CUBLAS_FILL_MODE_LOWER;
    transa = 'N';
    transb = 'C';
    uplo = 'L';
    start = min_n;
    end = max_n;
    alpha_d = alpha.real();
    beta_d = beta.real();


/*  Run a series of test with the increasing number of rows of C offloaded to the GPU  */

    printf("\n");
    printf("dev    GPU_N   CPU_N   K    Time(sec) GFLOPS\n");
    printf("%s\n", std::string(45,'-').c_str());
    best_time = 1.0e3;
    for( i = start; i <= end; i+=step){

        GPU_n = i;
        CPU_n = max_n - GPU_n;

        /* Compute GPU part */
        read_clock(&t_start);

        status = cublasXtZherkx(handle, cu_uplo, cu_trans, GPU_n, max_k, (const cuDoubleComplex*)&alpha, (const cuDoubleComplex *)A, lda, (const cuDoubleComplex *)B, ldb, &beta_d, (cuDoubleComplex *)C, ldc);
        if(status != CUBLAS_STATUS_SUCCESS){
            printf("Error in %s:%d CublasXtZherkx error! %d\n", __FILE__, __LINE__, status);
        }
        
        cudaThreadSynchronize();
        read_clock(&t_end);
        gpu_t = elapsed_time(&t_start, &t_end);
        gpu_gflops = (4.0 * (double)GPU_n * (double)GPU_n * (double)max_k)/(gpu_t*1e3);

        /* Compute CPU part */
        read_clock(&t_start);

        offset = GPU_n;
#ifdef MKL
        zgemmt_(&uplo, &transa, &transb, &CPU_n, &max_k, &alpha, &A[offset], &lda, &B[offset], &ldb, &beta, &C[GPU_n + GPU_n*ldc], &ldc);
#else
        zgemm_(&transa, &transb, &CPU_n, &CPU_n, &max_k, &alpha_d, &A[offset], &lda, &B[offset], &ldb, &beta_d, &C[GPU_n + GPU_n*ldc], &ldc);
#endif
        zgemm_(&transa, &transb, &CPU_n, &GPU_n, &max_k, &alpha_d, &A[offset], &lda, B, &ldb, &beta_d, &C[GPU_n], &ldc);

        read_clock(&t_end);
        cpu_t = elapsed_time(&t_start, &t_end);

    #ifdef MKL
        cpu_gflops = (4.0 * CPU_n * max_k * (CPU_n + 2.0 * GPU_n))/(cpu_t*1e3);
    #else
        cpu_gflops = (8.0 * CPU_n *max_k * (GPU_n + CPU_n))/(cpu_t*1e3);
    #endif

        printf("    %6d   %6d %6d\n", GPU_n, CPU_n, max_k);
        printf("GPU                        %7.4f   %7.4f\n", gpu_t*1e-6, gpu_gflops);
        printf("CPU                        %7.4f   %7.4f\n", cpu_t*1e-6, cpu_gflops);

        exec_time = max(gpu_t, cpu_t)*1e-6;
        if( exec_time < best_time ){
	        *ratio = (double)max_n / (double)GPU_n;
	        best_time = exec_time;
	    }
    }

    return;
}


void ratio_zher2k( cublasXtHandle_t handle, int min_n, int max_n, int max_k, int step, std::complex<double> *A, int lda, std::complex<double> *B, int ldb, std::complex<double> *C, int ldc, double *ratio )
{

    int    ndevices, devices[MAX_DEV], i;
    std::complex<double> alpha = 1.66 - 0.43i;
    std::complex<double> beta = -3.14 + 1.97i;
    double alpha_d, beta_d;
    char   transa, transb, uplo;
    int    n, k, GPU_n, CPU_n, offset;
    int    start, end;
    double gpu_gflops, cpu_gflops;
    double exec_time, best_time;

    /*  CUBLAS */
    cublasOperation_t cu_trans;
    cublasFillMode_t cu_uplo;
    cublasStatus_t status;
    
    /*  Timing routines  */
    struct timeval t_start, t_end;
    long int gpu_t, cpu_t;

    /*  Get number of cuda devices */
    cudaGetDeviceCount( &ndevices );

    /*  Set devices */
    for( i = 0; i < ndevices; i++)
        devices[i] = i;

    /*  Set parameters to zherk calls */
    cu_trans = CUBLAS_OP_N;
    cu_uplo = CUBLAS_FILL_MODE_LOWER;
    transa = 'N';
    transb = 'C';
    uplo = 'L';
    start = min_n;
    end = max_n;
    alpha_d = alpha.real();
    beta_d = beta.real();


/*  Run a series of test with the increasing number of rows of C offloaded to the GPU  */
 
    printf("\n");
    printf("dev    GPU_N   CPU_N   K    Time(sec) GFLOPS\n");
    printf("%s\n", std::string(45,'-').c_str());
    best_time = 1.0e3;

    for( i = start; i <= end; i+=step){

        GPU_n = i;
        CPU_n = max_n - GPU_n;

        /* Compute GPU part */
        read_clock(&t_start);

        status = cublasXtZher2k(handle, cu_uplo, cu_trans, GPU_n, max_k, (const cuDoubleComplex*)&alpha, (const cuDoubleComplex *)A, lda, (const cuDoubleComplex *)B, ldb, &beta_d, (cuDoubleComplex *)C, ldc);
        if(status != CUBLAS_STATUS_SUCCESS){
            printf("Error in %s:%d CublasXtZher2k error! %d\n", __FILE__, __LINE__, status);
        }
        
        cudaThreadSynchronize();
        read_clock(&t_end);
        gpu_t = elapsed_time(&t_start, &t_end);
        gpu_gflops = (8.0 * (double)GPU_n * (double)GPU_n * (double)max_k + 4.0*GPU_n)/(gpu_t*1e3);

        /* Compute CPU part */
        read_clock(&t_start);

        offset = GPU_n;

        zher2k_(&uplo, &transa, &CPU_n, &max_k, &alpha, &A[offset], &lda, &B[offset], &ldb, &beta_d, &C[GPU_n + GPU_n*ldc], &ldc);
        zgemm_(&transa, &transb, &CPU_n, &GPU_n, &max_k, &alpha_d, &B[offset], &ldb, A, &lda, &beta_d, &C[GPU_n], &ldc);
        zgemm_(&transa, &transb, &CPU_n, &GPU_n, &max_k, &alpha_d, &A[offset], &lda, B, &ldb, &beta_d, &C[GPU_n], &ldc);
	   
        read_clock(&t_end);
        cpu_t = elapsed_time(&t_start, &t_end);
        cpu_gflops = (8.0 * CPU_n * max_k * (CPU_n + 2.0*GPU_n) + 4.0*CPU_n)/(cpu_t*1e3);

        printf("    %6d   %6d %6d\n", GPU_n, CPU_n, max_k);
        printf("GPU                        %7.4f   %7.4f\n", gpu_t*1e-6, gpu_gflops);
        printf("CPU                        %7.4f   %7.4f\n", cpu_t*1e-6, cpu_gflops);

        exec_time = max(gpu_t, cpu_t)*1e-6;
        if( exec_time < best_time ){
            *ratio = (double)max_n / (double)GPU_n;
            best_time = exec_time;
        }
    }
    
    return;
}
