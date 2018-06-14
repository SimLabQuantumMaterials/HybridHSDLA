#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <omp.h>

#include "cx_double.h"
extern "C" {
#include "timing.h"
}
#include "types4stripped.h"
#include "blas.h"
#include "lapack.h"
#include "wrappers.h"

#include <mkl_service.h>
#ifdef GPU
    #include <cuda_runtime.h>
#endif

/* Reference from Markus Hrywniak
void apply_u_norm( int L, int a, int G, cx_double *B, int *lmaxs, int max_lmax, double *u_norms )
{
    //#pragma omp parallel for num_threads(16) schedule(dynamic,1)
    for (int at = 0; at < a; ++at) {
        int lmax = lmaxs[ at ];
        for (int l = 0; l <= lmax; ++l) {
            double sqrt_dE_u_norm = u_norms[ at*(max_lmax+1) + l ]; //sqrt(u_match(l, at).dE_u_norm);
            for (int m = -l; m <= l; ++m) {
                int lm = l*(l+1) + m;
                for ( int g = 0; g < G; g++ ) {
                    B[ lm + at*L + L*a*g ] *= sqrt_dE_u_norm;
                }
            }
        }
    }
}
*/

void apply_u_norm( int L, int a, int G, cx_double *B, int *lmaxs, int max_lmax, double *u_norms )
{
    int nths = omp_get_max_threads();
    double *linearized_norms = (double *) malloc (L*a*sizeof(double));
    // #pragma omp parallel for num_threads(nths)
    for (int at = 0; at < a; ++at) {
        int lmax = lmaxs[ at ];
        for (int l = 0; l <= lmax; ++l) {
            double sqrt_dE_u_norm = u_norms[ at*(max_lmax+1) + l ];
            int lm0 = l*(l+1) - l;
            int lme = lm0 + 2*l;
            for (int lm = lm0; lm <= lme; lm++) {
                linearized_norms[ lm + at*L ] = sqrt_dE_u_norm;
            }
        }
    }
    int lda = L*a;
    #pragma omp parallel for num_threads(nths)
    for( int col = 0; col < G; col++ )
        for( int row = 0; row < L*a; row++ )
            B[col*lda + row] *= linearized_norms[row];
    free( linearized_norms );
}

void apply_u_norm_s( int L, int a, int G, cx_double *B, cx_double *C, int *lmaxs, int max_lmax, double *u_norms )
{
    int nths = omp_get_max_threads();
    double *linearized_norms = (double *) malloc (L*a*sizeof(double));
    // #pragma omp parallel for num_threads(nths)
    for (int at = 0; at < a; ++at) {
        int lmax = lmaxs[ at ];
        for (int l = 0; l <= lmax; ++l) {
            double sqrt_dE_u_norm = u_norms[ at*(max_lmax+1) + l ];
            int lm0 = l*(l+1) - l;
            int lme = lm0 + 2*l;
            for (int lm = lm0; lm <= lme; lm++) {
                linearized_norms[ lm + at*L ] = sqrt_dE_u_norm;
            }
        }
    }
    int lda = L*a;
    #pragma omp parallel for num_threads(nths)
    for( int col = 0; col < G; col++ )
        for( int row = 0; row < L*a; row++ )
            C[col*lda + row] = B[col*lda + row] * linearized_norms[row];
    free( linearized_norms );
}

#define PTS printf("%5.2f ", omp_get_wtime() - fn_start_time);
#define TIC fn_tic = omp_get_wtime();
#define TOC(name, work, magnitude, unit) do { fn_tic = omp_get_wtime() - fn_tic; printf("%15s %5.2fs elapsed %8.2f %s\n", name, fn_tic, (work) / fn_tic / magnitude, unit); } while(0)

/*
 *                                  alive   needs   destroys
 * A, B = construct()
 * H = A* T B + B* T A + B* T B     A, B    A, B    A
 * A = (re)construct()              B
 * S = A* A + B* B                  A, B    A, B    B
 * H += A* TA                       A,      A       A, X
 */
void create_HS2_stripped_v2( int L, int a, int G,
                   cx_double *A, cx_double *B,
                   cx_double *H, cx_double *S,
                   myTMat *T,
                   int *l_maxs,
                   int max_lmax,
                   double *u_norms, 
                   std::string method,
                   cx_double *W )
{
    double fn_tic = 0;
    double fn_start_time = omp_get_wtime();
    // constants
    cx_double c0(0, 0);
    cx_double c1(1, 0);
    cx_double cp5(.5, 0);
    double dONE  = 1.0,
           dZERO = 0.0;

    // timings
    struct timeval t0, t1, t0op, t1op;
    long int elapsed_us;

    // Temporary matrices and pointers to matrices
    cx_double *X;
    cx_double *Tmp = (cx_double *) malloc ( L * G * sizeof(cx_double) );

    // leading dimensions
    int ldA = L * a;
    int ldB = L * a;
    int ldH = G;
    int ldS = G;
    int ldTmp = L;

    // Problem dimensions
    printf("[L] %d\n", L);
    printf("[A] %d\n", a);
    printf("[G] %d\n", G);


    // Allocate temporary array (workspace) for holding intermediate results from A or B
//    cx_double *W;
//#ifdef GPU
//    if( cudaHostAlloc( (void **)&W, (size_t)L * a * G * sizeof(cx_double), 0 ) != cudaSuccess )
//	printf("[ERROR] Couldn't allocate pinned memory for W\n");
//#else
//    W = (cx_double *) malloc ( L * a * G * sizeof(cx_double) );   
//#endif

    //TIC;
    read_clock(&t0);

    // S = A* A + B* B
    //
    // S = A* A
    TIC;
    int La = L * a;
    read_clock(&t0op);
    dev_zherk_("L", "C", &G, &La, &dONE, A, &ldA, &dZERO, S, &ldS);
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[S: zherk 1] %.3f secs\n", elapsed_us*1e-6);
    TOC("S zherk 1", 4.0 * G * 1.0 * G * 1.0 * La, 1.0e9, "GFLOPS/s");

    // apply u norm to B
    read_clock(&t0op);
//    apply_u_norm(L, a, G, B, l_maxs, max_lmax, u_norms );
    apply_u_norm_s(L, a, G, B, W, l_maxs, max_lmax, u_norms );
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[S: unorm] %.3f secs\n", elapsed_us*1e-6);
    printf("    unorm   %4.2fs elapsed  %.2f GFLOPS/s\n", elapsed_us*1e-6,
                                                       2 * L * a * G * 1e-9 / (elapsed_us*1e-6)); // rougly

    // S += B* B
    TIC;
    read_clock(&t0op);
    dev_zherk_("L", "C", &G, &La, &dONE, W, &ldB, &dONE, S, &ldS);
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[S: zherk 2] %.3f secs\n", elapsed_us*1e-6);
    TOC("S zherk 2", 4.0 * G * 1.0 * G * 1.0 * La, 1.0e9, "GFLOPS/s");


    // X stored in W
    X = W;
    int ldX = ldA;
    int X_height = 0;


    read_clock(&t0op);

#if 1
    for (size_t atom = 0; atom < a; atom++) {
       // Tmp = 0
       bzero( Tmp, L * G * sizeof(cx_double) );

       // Tmp = T_BB[atom] B[atom]
       cx_double *T_BB_atom = const_cast<cx_double*>(T[atom].bb);
       int T_BB_atom_size = T[atom].bb_size;
       int ldT_BB_atom = T_BB_atom_size;
       cx_double *B_atom = B + atom * L;
       zhemm_("L", "L", (int*)&T_BB_atom_size, &G, &cp5, T_BB_atom, &ldT_BB_atom, B_atom, &ldB, &c0, Tmp, &ldTmp);

       // Tmp += T_AB[atom] A[atom]
       cx_double *T_AB_atom = const_cast<cx_double*>(T[atom].ab);
       int T_AB_atom_size = T[atom].ab_size;
       int ldT_AB_atom = T_AB_atom_size;
       cx_double *A_atom = A + atom * L;
       zgemm_("C", "N", &T_AB_atom_size, &G, &T_AB_atom_size, &c1, T_AB_atom, &ldT_AB_atom, A_atom, &ldA, &c1, Tmp, &ldTmp);

       // [ALGORITHMIC] Notice that T_AB_atom_size, T_BB_atom_size might be smaller than L (size of Tmp)
       // That is why we zero it out and then below take the max for next step (to copy only
       // the nonzero parts of Tmp and avoiding multiplications by zero later on)

       // X[X_height] = Tmp
       int Tmp_size = (T_AB_atom_size > T_BB_atom_size) ? T_AB_atom_size : T_BB_atom_size;
       cx_double *X_X_height = X + X_height;
       zlacpy_("A", &Tmp_size, &G, Tmp, &ldTmp, X_X_height, &ldX);

       if (X_height != atom * L) {
          // Tmp = B_atom
          zlacpy_("A", &Tmp_size, &G, B_atom, &ldB, Tmp, &ldTmp);

          // B[X_height] = Tmp
          cx_double *B_X_height = B + X_height;
          //zlacpy_("A", &Tmp_size, &G, Tmp, &ldTmp, B_X_height, &ldTmp);
          zlacpy_("A", &Tmp_size, &G, Tmp, &ldTmp, B_X_height, &ldB);
       }

       X_height += Tmp_size;
    }
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[Atoms 1] %.3f secs\n", elapsed_us*1e-6);
    printf("    Atoms 1   %4.2fs elapsed  %.2f GFLOPS/s\n", elapsed_us*1e-6,
                                                             8.0 * 2 * L * L * a * G * 1e-9 / (elapsed_us*1e-6)); // rougly
#endif

#if 0
        //read_clock(&t0op);
        mkl_set_num_threads(1);
#pragma omp parallel num_threads(NTHS)
    {
        cx_double *Tmp = (cx_double *) malloc ( L * G * sizeof(cx_double) );
        int ldTmp = L;

#pragma omp for
        for (size_t atom = 0; atom < a; atom++) {
            // Tmp = 0
            bzero( Tmp, L * G * sizeof(cx_double) );

            // Tmp = T_BB[atom] B[atom]
            cx_double *T_BB_atom = const_cast<cx_double*>(T[atom].bb);
            int T_BB_atom_size = T[atom].bb_size;
            int ldT_BB_atom = T_BB_atom_size;
            cx_double *B_atom = B + atom * L;
            zhemm_("L", "L", (int*)&T_BB_atom_size, &G, &cp5, T_BB_atom, &ldT_BB_atom, B_atom, &ldB, &c0, Tmp, &ldTmp);

            // Tmp += T_AB[atom] A[atom]
            cx_double *T_AB_atom = const_cast<cx_double*>(T[atom].ab);
            int T_AB_atom_size = T[atom].ab_size;
            int ldT_AB_atom = T_AB_atom_size;
            cx_double *A_atom = A + atom * L;
            zgemm_("C", "N", &T_AB_atom_size, &G, &T_AB_atom_size, &c1, T_AB_atom, &ldT_AB_atom, A_atom, &ldA, &c1, Tmp, &ldTmp);

            // [ALGORITHMIC] Notice that T_AB_atom_size, T_BB_atom_size might be smaller than L (size of Tmp)
            // That is why we zero it out and then below take the max for next step (to copy only
            // the nonzero parts of Tmp and avoiding multiplications by zero later on)

            // X[X_height] = Tmp
            int Tmp_size = (T_AB_atom_size > T_BB_atom_size) ? T_AB_atom_size : T_BB_atom_size;
            cx_double *X_X_height = X + X_height;
            zlacpy_("A", &Tmp_size, &G, Tmp, &ldTmp, X_X_height, &ldX);

            if (X_height != atom * L) {
                // Tmp = B_atom
                zlacpy_("A", &Tmp_size, &G, B_atom, &ldB, Tmp, &ldTmp);

                // B[X_height] = Tmp
                cx_double *B_X_height = B + X_height;
                zlacpy_("A", &Tmp_size, &G, Tmp, &ldTmp, B_X_height, &ldTmp);
            }

            X_height += Tmp_size;
        }
        free( Tmp );
    }
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[Atoms 1] %.2f secs\n", elapsed_us*1e-6);
        mkl_set_num_threads(NTHS);

        X_height = L*a;
#endif
        // H = B* X + X* B
    printf("[X_height] %d\n", X_height);
    TIC;
    read_clock(&t0op);
    double zero = 0.0;
    dev_zher2k_("L", "C", &G, &X_height, &c1, B, &ldB, X, &ldX, &zero, H, &ldH);
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[H1: zher2k] %.3f secs\n", elapsed_us*1e-6);
    TOC("H zher2k", 8 * G * 1.0 * G * 1.0 * X_height, 1.0e9, "GFLOPS/s");


    //TIC;
    // H += A* T_AA A
    read_clock(&t0op);
    {
        // X stored in B
        X = B;
        int ldX = ldB;

        int nhpd = 0;
        int X_nhpd_height = 0;
        int X_hpd_height = 0;

        // [ALGORITHMIC] Notice that below the nonhpd part overwrites from top to bottom
        // matrix X, while the hpd part does from bottom to top
        for (size_t atom = 0; atom < a; atom++) {
            // Cholbuff = Chol(T_AA[atom])
            cx_double *T_AA_atom = const_cast<cx_double*>(T[atom].aa);
            int T_AA_atom_size = T[atom].aa_size;
            int ldT_AA_atom = T_AA_atom_size;
            int info;

            cx_double *A_atom = A + atom * L;
            int idx = X_nhpd_height;

            // A[idx] = A[atom]
            cx_double *A_idx = A + idx;
            if (A_idx != A_atom)
                //zlacpy_("A", &T_AA_atom_size, &G, A_atom, &ldA, A_idx, &ldA);
            {
                // Tmp = A[atom]
                zlacpy_("A", &T_AA_atom_size, &G, A_atom, &ldA, Tmp, &ldTmp);
                // A[idx] = Tmp
                cx_double *A_height = A + X_nhpd_height;
                zlacpy_("A", &T_AA_atom_size, &G, Tmp, &ldTmp, A_height, &ldA);
            }

            // X[idx] = T_AA[atom] A[idx]
            cx_double *X_idx = X + idx;
            zhemm_("L", "L", &T_AA_atom_size, &G, &c1, T_AA_atom, &ldT_AA_atom, A_idx, &ldA, &c0, X_idx, &ldA);

            X_nhpd_height += T_AA_atom_size;
            //}
        }
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[Atoms 2] %.3f secs\n", elapsed_us*1e-6);
        printf("    Atoms 2   %4.2fs elapsed  %.2f GFLOPS/s\n", elapsed_us*1e-6,
                                                             8.0 * L * L * a * G * 1e-9 / (elapsed_us*1e-6)); // rougly

        // H += A[nhpd]* X[nhpd]
        printf("[X_nhpd_height] %d\n", X_nhpd_height);
        TIC;
        read_clock(&t0op);
        cx_double *A_nhpd = A;
        cx_double *X_nhpd = X;
        dev_zherkx_("L", "C", &G, &X_nhpd_height, &c1, A_nhpd, &ldA, X_nhpd, &ldX, &dONE, H, &ldH);
        // Just to check that automatic offloading for phi works:
        //zgemm_("T", "N", &G, &G, &X_nhpd_height, &c1, A_nhpd, &ldA, X_nhpd, &ldX, &c1, H, &ldH);
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[H2 zherkx] %.3f secs\n", elapsed_us*1e-6);
        TOC("H zherkx", 4 * G * 1.0 * G * X_nhpd_height, 1.0e9, "GFLOPS/s");
    }
    //TOC("H += A* T_AA A                             ");
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[Compute time] %.3f secs\n", elapsed_us*1e-6);

    free(Tmp);
//#ifdef GPU
//    cudaFreeHost(W);
//#else
//    free(W);
//#endif

}

