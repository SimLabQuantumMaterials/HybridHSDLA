/*
* The program computes how many rows/columns of matrix C in BLAS-3 operations have to be offloaded to the GPU to achieve the maximum performance when load-balancing workload between GPUs and CPUs.
* 
* The program returns the ratio of the rows of C offloaded to the GPUs agains the total number of rows (i.e. ratio is computed as n/number_of_gpu_rows). The number of rows to be offloaded to the GPU
* is then computed for each case (with row dimension, e.g. 'm') as number_of_gpu_rows = m/ratio.
*
* The ratio is computed for each BLAS-3 kernel. Currently, only zherk, zher2k and zherkx are implemented.
*
*/

#include <stdio.h>
#include <string.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cublasXt.h>
#include <cuda.h>
#include <complex>
#include <omp.h>

#ifdef MKL
  #include <mkl_service.h>
#else
  #include <cblas.h>
#endif

/* Timing routines */
#include "timings.h"
#include "kernels.h"

#define MIN_DIM 1000
#define MAX_DIM 5000
#define MAX_K   10000
#define STEP_BLK 128
#define MAX_DEV 16


/*  External functions  */
extern "C"{
	void zlarnv_(int*, int*, int*, std::complex<double>*);
}


/* Auxiliary function  */
void print_usage(char *progname){
	printf("The program computes how many rows/columns of matrix C in BLAS-3 operations have to be offloaded to the GPU to achieve the maximum performance when load-balancing workload between GPUs and CPUs.\n" 
		"The return value is a ratio that determines the percentage of rows of C to be offloaded to the GPUs (the gpu part is computed as n/ratio). The ratio is computed for each BLAS-3 kernel.\n" 
		"Currently, only zherk, zher2k and zherkx are implemented.\n");
	printf("\nUsage: ");
	printf("%s -n min max -k k -s step -g gpus -t thr -h \n\n", progname);
	printf(" Parameters:\n");
	printf(" -n      min and max number of rows and columns of matrix C and rows/columns of op(A)/op(B) (default = min: %d max: %d)\n", MIN_DIM, MAX_DIM);
	printf(" -k      number of columns of matrix op(A) and rown op(B) (default = %d)\n", MAX_K);
	printf(" -s      step size for iteratively increasing size of matrix offloaded to the GPU (default = %d)\n", STEP_BLK);
	printf(" -g      number of GPUs (default = all gpus)\n");
	printf(" -t      number of threads for a hybrid variant (default = all cpus)\n");
	printf(" -h      prints this help message)\n");
	printf("\n");
}

int main(int argc, char **argv)
{
	std::complex<double> *A = NULL, *B = NULL, *C = NULL;
	int 	lda, ldb, ldc;
	int 	max_n, min_n, k, step, nthrs;
	int	i;
	double 	zherk_r, zher2k_r, zherkx_r;
	int 	nelem, idist = 2;
	int 	iseed1[] = {1,11,7,1}, iseed2[] = {3,7,13,13}, iseed3[] = {3,11,13,1};

	/*  CUBLAS parameters */
	cublasStatus_t stat;
	cudaError_t c_stat;
	cublasXtHandle_t handle;
	int ndevs, devices[MAX_DEV];


	/* Set default values */
	min_n = MIN_DIM;
	max_n = MAX_DIM;
	k = MAX_K;
	step = STEP_BLK;
	nthrs = 0;
	ndevs = 0;

	/* Read input parameters */
	for( i = 1; i < argc; i++ ){
		if(strcmp(argv[i],"-n") == 0 && i+2 < argc){
			min_n = atoi(argv[++i]);
			max_n = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-k") == 0 && i+1 < argc){
			k = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-t") == 0 && i+1 < argc){
			nthrs = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-s") == 0 && i+1 < argc){
			step = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-g") == 0 && i+1 < argc){
			ndevs = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-h") == 0){
			print_usage(argv[0]);
			exit(0);
		}else{
			printf("Unknown parameter: %s\n", argv[i]);
			print_usage(argv[0]);
			exit(-1);
		}
	}

	/* Initialize CublasXt */
	stat = cublasXtCreate(&handle);
	if (stat != CUBLAS_STATUS_SUCCESS){
		fprintf(stderr, "[%s:%d] CUBLASXT initialization error\n", __FILE__, __LINE__);
		return -2;
    	}

	/*  If number of threads is not set, then use all available threads  */
	if(nthrs == 0){
		nthrs = omp_get_max_threads();
	}

	/*  If number of GPU devices is not set, the use all devices of the system  */
	if(ndevs == 0){
		c_stat = cudaGetDeviceCount(&ndevs);
		if( stat != CUBLAS_STATUS_SUCCESS){
			fprintf(stderr, "[%s:%d] Error getting number of devices\n", __FILE__, __LINE__);
			fprintf(stderr, "Error code: %d\n", c_stat);
			return -2;
		}
	}

	/*  Print test configuration  */
        printf("%s\n", std::string(45,'=').c_str());
	printf("CONFIGURATION\n");
	printf("\n");
	printf("N range: %d - %d\n", min_n, max_n);	
	printf("K:       %d\n", k);	
	printf("Step:    %d\n", step);	
	printf("Threads: %d\n", nthrs);	
	printf("GPUs:    %d\n", ndevs);	
        printf("%s\n", std::string(45,'=').c_str());

	/* Set leading dimensions of the testing matrices */
    	lda = ldb = ldc = max_n;

	/*  Allocate matrices  */
	if ( cudaMallocHost( (void**)&A, lda*k*sizeof(std::complex<double>) ) != cudaSuccess ){
    		printf("[%s: %d]: Error in allocating pinned memory for A!\n", argv[0], __LINE__);
    	}
	if ( cudaMallocHost( (void**)&B, ldb*k*sizeof(std::complex<double>) ) != cudaSuccess ){
		printf("[%s: %d]: Error in allocating pinned memory for B\n", argv[0], __LINE__);
	}
	if ( cudaMallocHost( (void**)&C, ldc*max_n*sizeof(std::complex<double>) ) != cudaSuccess ){
		printf("[%s: %d]: Error in allocating pinned memory for C\n", argv[0], __LINE__);
	}

	/*  Generate matrices with the random numbers  */
        nelem = max_n*k;
	zlarnv_(&idist, iseed1, &nelem, A);
	zlarnv_(&idist, iseed2, &nelem, B);
        nelem = max_n*max_n;
	zlarnv_(&idist, iseed3, &nelem, C);

	/*  Set GPU devices  */
	for( int i = 0; i < ndevs; i++)
		devices[i] = i;

	/*  Set number of GPU devices for CublasXt  */
	stat = cublasXtDeviceSelect(handle, ndevs, devices);

	/*  Set number of threads  */
	#ifdef MKL
    	    mkl_set_num_threads(nthrs);
	#else
    	    openblas_set_num_threads(nthrs);
	#endif

	/* Compute ratio for zherk */
	printf("\nTesting ZHERK...\n");
	ratio_zherk(handle, min_n, max_n, k, step, A, lda, C, ldc, &zherk_r);
	printf("Ratio for zherk is: %.4f\n", zherk_r);

	/* Compute ratio for zher2k */
	printf("\n Testing ZHER2K...\n");
	ratio_zher2k(handle, min_n, max_n, k, step, A, lda, B, ldb, C, ldc, &zher2k_r);
	printf("Ratio for zher2k is: %.4f\n", zher2k_r);

	/* Compute ratio for zherkx */
	printf("\n Testing ZHERKX...\n");
	ratio_zherkx(handle, min_n, max_n, k, step, A, lda, B, ldb, C, ldc, &zherkx_r);
	printf("Ratio for zherkx is: %.4f\n", zherkx_r);

		
	/* Close CublasXt */
	stat = cublasXtDestroy(handle);
        if (stat != CUBLAS_STATUS_SUCCESS){
        	fprintf(stderr, "[%s] CublasXtDestroy error\n", argv[0]);
 	      	return -1;
    	}

	printf("\n");

    /*  Deallocate arrays on the host  */
	cudaFreeHost(A);
	cudaFreeHost(B);
	cudaFreeHost(C);

	return EXIT_SUCCESS;
}
