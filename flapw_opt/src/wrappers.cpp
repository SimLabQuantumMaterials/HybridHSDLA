#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <errno.h>

#include "cx_double.h"
extern "C" {
    #include "timing.h"
}
#include "blas.h"
#include "types4stripped.h"
#include "wrappers.h"

#ifdef GPU
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    #include <cublasXt.h>
    #include <cuda_profiler_api.h>
#endif
#include <omp.h>
#include <mkl_service.h>
#if defined HS_CUDA ||  HS_HSTREAMS
    #include "hetero-mic-gpu/source.h"
#endif

#ifdef GPU
static cublasXtHandle_t handle;
static int *devices;
#endif

#ifdef OOC
    static cublasHandle_t handleCublas;
    int block_size;
    double ratio[3];
    int nthreads;

extern "C"{
    int zherkx_hyb( cublasXtHandle_t handle, char uplo, char trans, int n, int k, std::complex<double> alpha, std::complex<double> *A, int lda, std::complex<double> *B, int ldb, std::complex<double> beta, std::complex<double> *C, int ldc, int nb, double ratio, int nthreads );
    int zherk_hyb( cublasXtHandle_t handle, char uplo, char trans, int n, int k, double alpha, std::complex<double> *A, int lda, double beta, std::complex<double> *C, int ldc, int nb, double ratio, int nthreads );
    int zher2k_hyb( cublasXtHandle_t handle, char uplo, char trans, int n, int k, std::complex<double> alpha, std::complex<double> *A, int lda, std::complex<double> *B, int ldb, double beta, std::complex<double> *C, int ldc, int nb, double ratio, int nthreads );
    void get_ratio( int nthreads, double *ratio );
    void LOG_CUBLAS_STATUS(cublasStatus_t status, char *);
}
#endif

int get_env_int( char *var )
{
    char *value_str, *end;
    long int value;

    value_str = getenv( var );
    if ( value_str == NULL )
    {
        printf("%s unset\n", var);
        value = -1;
    }
    else
    {
        errno = 0;
        value = strtol( value_str, &end, 10 );
        if ( errno == ERANGE || value_str == end && strcmp( value_str, "" ) )
        {
            fprintf( stderr, "[WARNING] Incorrect value for $%s: '%s'\n", var, value_str );
            value = -1;
        }
    }
    return (int)value;
}

void get_ratio( int dev_count, double* ratio)
{
    char var[] = "HTCON_RATIO_FILE";
    char *ratio_file;
    FILE *stream;
    char *line_buffer;
    size_t len = 0;
    ssize_t read;
    int counter = 0;

    ratio_file = getenv( var );
    if( ratio_file == NULL) 
	printf("%s unset\n", var);

    stream = fopen(ratio_file, "r");
    if ( stream == NULL ){
	fprintf(stderr, "Error opening file: %s (%s:%d)\n", ratio_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
    }

    /* Check if the ratio file has the number of lines >= dev_count */
    
    // TODO

    line_buffer = (char *)malloc(100 * sizeof(char));

    /* Read dev_count-th line from the file and save ratio values */
    while(read = getline(&line_buffer, &len, stream) != -1){
	if(++counter == dev_count){
	    sscanf(line_buffer, "%lf %lf %lf", &ratio[0], &ratio[1], &ratio[2]);
	    break;
        }
    }
    fclose(stream);
    free(line_buffer);

}

int dev_init( )
{
#if defined HS_CUDA || defined HS_HSTREAMS
    hs_init(NULL, NULL, 0, NULL, NULL, 0);
#elif defined GPU
    cublasStatus_t status;
    int i, dev_count;

    dev_count = get_env_int("HTCON_NDEVS");
    if ( dev_count == -1 )
        cudaGetDeviceCount( &dev_count );
    printf("[DEVICES] %d devices\n", dev_count );
    devices = (int *) malloc( dev_count * sizeof(int) );
    for( i = 0; i < dev_count; i++ )
        devices[i] = i;

    status = cublasXtCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        fprintf(stderr, "[myCUDAwrapper] CUBLASXT initialization error\n");
        return -2;
    }

    LOG_CUBLAS_STATUS(cublasXtDeviceSelect(handle, dev_count, devices), "cublasXtDevices");
#ifdef OOC
    nthreads = omp_get_max_threads();

    get_ratio(dev_count, ratio);

    printf("Ratios: [%f, %f, %f]\n", ratio[0], ratio[1], ratio[2]);

#endif
    int xt_tile_size = get_env_int("HTCON_XT_TS");
    if (xt_tile_size != -1)
        cublasXtSetBlockDim( handle, xt_tile_size);
#endif
    return 0;
}

int dev_finalize( )
{
#if defined HS_CUDA || defined HS_HSTREAMS
    hs_fini();
#elif defined GPU
    cublasStatus_t status;
        
    status = cublasXtDestroy(handle);
    cudaDeviceReset();
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        fprintf(stderr, "[myCUDAwrapper] Shutdown error\n");
        return -1;
    }

    free( devices );
#endif

    return 0;
}

void dev_zherkx_( char *uplo, char *transa, int *n, int *k,
                  cx_double *alpha, cx_double *A, int *lda,
                                    cx_double *B, int *ldb, 
                     double *beta,  cx_double *C, int *ldc )
{
    int m[1];
    m[0] = *n;
    char transb[1];
    if (transa[0] == 'N')
        transb[0] = 'C';
    else
        transb[0] = 'N';
#if defined HS_CUDA || defined HS_HSTREAMS
    hs_zherkx(uplo, transa, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#elif defined GPU
    printf("[INFO] Offloading zherkx to gpu\n");
    cublasFillMode_t cu_uplo   = uplo[0] == 'L' ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
    cublasOperation_t cu_transa = transa[0] == 'N' ? CUBLAS_OP_N : 
                                  transa[0] == 'T' ? CUBLAS_OP_T : CUBLAS_OP_C;
#ifdef OOC
    int nthreads = omp_get_max_threads();
    std::complex<double> beta_complex = *beta;
    printf("CALLING HERKX_HYB\n");
    zherkx_hyb( handle, uplo[0], transa[0], *n, *k, *alpha, A, *lda, B, *ldb, beta_complex, C, *ldc, 0, ratio[2], nthreads );
#else
    printf("CALLING HERKX\n");
    cublasXtZherkx( handle, cu_uplo, cu_transa, *n, *k, 
            (cuDoubleComplex *)alpha, (cuDoubleComplex *)A, *lda,
            (cuDoubleComplex *)B, *ldb,
            beta,  (cuDoubleComplex *)C, *ldc );
#endif    
#else
    cx_double cx_beta(*beta, 0);
    zgemmt_(uplo, transa, transb, n, k, alpha, A, lda, B, ldb, &cx_beta, C, ldc);
#endif    
}

void dev_zherk_( char *uplo, char *trans,  int *n,  int *k,
              double *alpha, cx_double *A,  int *lda,
              double *beta,  cx_double *C,  int *ldc )
{
#if defined HS_CUDA || defined HS_HSTREAMS
    hs_zherk(uplo, trans, n, k, alpha, (cx_double*) A, lda, beta, (cx_double*) C, ldc);
#elif defined GPU
    printf("[INFO] Offloading zherk to gpu\n");
    cublasFillMode_t cu_uplo = uplo[0] == 'L' ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
    cublasOperation_t cu_trans = trans[0] == 'N' ? CUBLAS_OP_N : 
                                 trans[0] == 'T' ? CUBLAS_OP_T : CUBLAS_OP_C;
#ifdef OOC
    int nthreads = omp_get_max_threads();
    printf("CALLING ZHERK_HYB\n");
    zherk_hyb( handle, uplo[0], trans[0], *n, *k, *alpha, A, *lda, *beta, C, *ldc, 0, ratio[0], nthreads );
#else
    cublasXtZherk( handle, cu_uplo, cu_trans, *n, *k, 
                   alpha, (cuDoubleComplex *)A, *lda,
                   beta,  (cuDoubleComplex *)C, *ldc );
#endif
#else
    zherk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
#endif
}

void dev_zher2k_( char *uplo, char *trans,  int *n,  int *k,
                  cx_double *alpha,  cx_double *A,  int *lda,
                  cx_double *B,  int *ldb, 
                  double *beta, cx_double *C,  int *ldc )
{
#if defined HS_CUDA || defined HS_HSTREAMS
    hs_zher2k(uplo, trans, n, k, alpha, A, lda, B, ldb, (const double *)beta, C, ldc);
#elif defined GPU
    printf("[INFO] Offloading zher2k to gpu\n");
    cublasFillMode_t cu_uplo = uplo[0] == 'L' ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
    cublasOperation_t cu_trans = trans[0] == 'N' ? CUBLAS_OP_N : 
                                 trans[0] == 'T' ? CUBLAS_OP_T : CUBLAS_OP_C;
#ifdef OOC
    printf("CALLING ZHER2K_HYB\n");
    zher2k_hyb( handle, uplo[0], trans[0], *n, *k, *alpha, A, *lda, B, *ldb, *beta, C, *ldc, 0, ratio[1], nthreads );
#else
    cublasXtZher2k( handle, cu_uplo, cu_trans, *n, *k, 
                    (cuDoubleComplex *)alpha, (cuDoubleComplex *)A, *lda,
                                              (cuDoubleComplex *)B, *ldb,
                                        beta, (cuDoubleComplex *)C, *ldc );
#endif
#else
    zher2k_(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}
