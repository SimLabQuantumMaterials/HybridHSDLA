#ifndef _WRAPPERS
#define _WRAPPERS

#include "cx_double.h"

int dev_init( );
int dev_finalize( );

void dev_zgemm_( char *transa, char *transb, int *m, int *n, int *k,
                 cx_double *alpha, cx_double *A, int *lda,
                                   cx_double *B, int *ldb, 
                 cx_double *beta,  cx_double *C, int *ldc );

void dev_zhemm_( char *side, char *uplo, int *m, int *n,
                 cx_double *alpha, cx_double *A, int *lda,
                                   cx_double *B, int *ldb, 
                 cx_double *beta,  cx_double *C, int *ldc );

void dev_zherk_( char *uplo, char *trans, int *n, int *k,
                 double *alpha, cx_double *A, int *lda,
                 double *beta,  cx_double *C, int *ldc );

void dev_zher2k_( char *uplo, char *trans, int *n, int *k,
                  cx_double *alpha, cx_double *A, int *lda,
                                    cx_double *B, int *ldb, 
                     double *beta,  cx_double *C, int *ldc );

void dev_zherkx_( char *uplo, char *transa, int *n, int *k,
                  cx_double *alpha, cx_double *A, int *lda,
                                    cx_double *B, int *ldb, 
                     double *beta,  cx_double *C, int *ldc );

#endif // _WRAPPERS
