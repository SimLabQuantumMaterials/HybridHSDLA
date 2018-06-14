#ifndef _LAPACK
#define _LAPACK

#include "cx_double.h"

// BLAS/LAPACK signatures
extern "C" {
    void zlacpy_( const char *, const int *, const int *, const cx_double *, const int *, cx_double *, const int *);
    /*void zpotrf_( const char *, const int *, cx_double *, const int *, int *);*/
    void dlarnv_( const int *, int *, const int *, double *);
}

#endif // _LAPACK
