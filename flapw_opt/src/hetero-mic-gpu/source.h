#ifndef _SOURCE
#define _SOURCE

#include "cx_double.h"

// BLAS/LAPACK signatures
extern "C" {
    void hs_zgemm(const char *, const char *,  
                  const int *, const int *, const int *,
                  const cx_double *, const cx_double *, const int *,  
                                     const cx_double *, const int *,  
                  const cx_double *,       cx_double *, const int *);
    void hs_zhemm(const char *, const char *,  
                  const int *, const int *,  
                  const cx_double *, const cx_double *, const int *,  
                                     const cx_double *, const int *,  
                  const cx_double *,       cx_double *, const int *);
    void hs_zherk(const char *, const char *,  
                  const int *, const int *,  
                  const double *, const cx_double *, const int *,  
                  const double *,       cx_double *, const int *);
    void hs_zher2k(const char *, const char *,  
                   const int *, const int *,
                   const cx_double *, const cx_double *, const int *,  
                                      const cx_double *, const int *,  
                   const    double *,       cx_double *, const int *);
    void hs_zherkx(const char *, const char *,  
                   const int *, const int *,
                   const cx_double *, const cx_double *, const int *,  
                                      const cx_double *, const int *,  
                   const    double *,       cx_double *, const int *);
}

extern "C" {
    void hs_init(void * A, void * B, uint64_t szAB, void * H, void * S, uint64_t szHS);
    void hs_fini();
}

#endif // _SOURCE
