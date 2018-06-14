#ifndef _BLAS
#define _BLAS

#include "cx_double.h"

// BLAS/LAPACK signatures
extern "C" {
    void dgemm_(const char *, const char *,  
                const int *, const int *, const int *,
                const double *, const double *, const int *,  
                                const double *, const int *,  
                const double *,       double *, const int *);
    void zgemm_(const char *, const char *,  
                const int *, const int *, const int *,
                const cx_double *, const cx_double *, const int *,  
                                   const cx_double *, const int *,  
                const cx_double *,       cx_double *, const int *);
    void zhemm_(const char *, const char *,  
                const int *, const int *,  
                const cx_double *, const cx_double *, const int *,  
                                   const cx_double *, const int *,  
                const cx_double *,       cx_double *, const int *);
    void zherk_(const char *, const char *,  
                const int *, const int *,  
                const double *, const cx_double *, const int *,  
                const double *,       cx_double *, const int *);
    void zher2k_(const char *, const char *,  
                 const int *, const int *,
                 const cx_double *, const cx_double *, const int *,  
                                    const cx_double *, const int *,  
                 const    double *,       cx_double *, const int *);
    void ztrmm_(char *, char *, char *, char *, int *,  int *,  cx_double *,  cx_double *,  int *, cx_double *,  int *);
    void zgemmt_(const char *, const char *, const char *,
                 const int *, const int *,
                 const cx_double *, const cx_double *, const int *, 
                                    const cx_double *, const int *,
                 const cx_double *,       cx_double *, const int *);
    /*void zlacpy_( const char *, const int *, const int *, const cx_double *, const int *, cx_double *, const int *);*/
    /*void zpotrf_( const char *, const int *, cx_double *, const int *, int *);*/
}

#endif // _BLAS
