#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <complex>
#include <omp.h>
#include "types4stripped.h"

#include "myCUDAwrapper.h"
extern "C" {
    #include "timing.h"
}

#include <mkl_service.h>

#define NTHS 24

typedef std::complex<double> cx_double;

// declare some BLAS/LAPACK signatures
extern "C" {
    void zhemm_( char *,  char *,  int *,  int *,  cx_double *,  cx_double *,  int *,  cx_double *,  int *,  cx_double *, cx_double *,  int *);
    void zher2k_( char *,  char *,  int *,  int *,  cx_double *,  cx_double *,  int *,  cx_double *,  int *,  cx_double *, cx_double *,  int *);
    void zlacpy_( char *,  int *,  int *,  cx_double *,  int *, cx_double *,  int *);
    void zpotrf_( char *,  int *, cx_double *,  int *, int *);
    void ztrmm_( char *,  char *,  char *,  char *,  int *,  int *,  cx_double *,  cx_double *,  int *, cx_double *,  int *);
    void zherk_( char *,  char *,  int *,  int *,  double *,  cx_double *,  int *,  double *, cx_double *,  int *);
    void zgemm_( char *,  char *,  int *,  int *,  int *,  cx_double *,  cx_double *,  int *,  cx_double *,  int *,  cx_double *, cx_double *,  int *);
}

/*
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

void apply_u_norm( int L, int a, int G, cx_double *B, int *lmaxs, int max_lmax, double *u_norms )
{
#pragma omp parallel for num_threads(16) schedule(dynamic,1)
    for (int at = 0; at < a; ++at) {
        int lmax = lmaxs[ at ];
        for (int l = 0; l <= lmax; ++l) {
            double sqrt_dE_u_norm = u_norms[ at*(max_lmax+1) + l ]; //sqrt(u_match(l, at).dE_u_norm);
            int lm0 = l*(l+1);
            int lme = lm0 + 2*l;
            for ( int g = 0; g < G; g++ ) {
                for (int lm = lm0; lm <= lme; lm++) {
                    B[ lm + at*L + L*a*g ] *= sqrt_dE_u_norm;
                }
            }
        }
    }
}
*/

void apply_u_norm( int L, int a, int G, cx_double *B, int *lmaxs, int max_lmax, double *u_norms )
{
    double *linearized_norms = (double *) malloc (L*a*sizeof(double));
    for (int at = 0; at < a; ++at) {
        int lmax = lmaxs[ at ];
        for (int l = 0; l <= lmax; ++l) {
            double sqrt_dE_u_norm = u_norms[ at*(max_lmax+1) + l ]; //sqrt(u_match(l, at).dE_u_norm);
            int lm0 = l*(l+1) - l;
            int lme = lm0 + 2*l;
            for (int lm = lm0; lm <= lme; lm++) {
                //printf("Idx: %d", lm + at*L);
                //fflush(stdout);
                linearized_norms[ lm + at*L ] = sqrt_dE_u_norm;
            }
        }
    }
    int lda = L*a;
    #pragma omp parallel for num_threads(NTHS)
    for( int col = 0; col < G; col++ )
        for( int row = 0; row < L*a; row++ )
            B[col*lda + row] *= linearized_norms[row];
    free( linearized_norms );
}


/*
 *                                  alive   needs   destroys
 * A, B = construct()
 * H = A* T B + B* T A + B* T B     A, B    A, B    A
 * A = (re)construct()              B
 * S = A* A + B* B                  A, B    A, B    B
 * H += A* TA                       A,      A       A, X
 */
void create_HS2_stripped( int L, int a, int G,
                   cx_double *A, cx_double *B,
                   cx_double *H, cx_double *S,
                   myTMat *T,
                   int *l_maxs,
                   int max_lmax,
                   double *u_norms, 
                   std::string method, int tile_size )
{
    double fn_tic = 0;
#   define TIC fn_tic = omp_get_wtime();
#   define TOC(name, work, magnitude, unit) do { fn_tic = omp_get_wtime() - fn_tic; printf("%15s %5.2fs elapsed %8.2f %s\n", name, fn_tic, (work) / fn_tic / magnitude, unit); } while(0)
    // constants
    cx_double c0(0, 0);
    cx_double c1(1, 0);
    cx_double cp5(.5, 0);
    cx_double *X;

    double dONE  = 1.0,
           dZERO = 0.0;

    // timings
    struct timeval t0, t1, t0op, t1op;
    long int elapsed_us;

    //#define TIC tic = Clock::now()
    //#define TOC(msg) dt = diff_as_msecs(tic); LOG(INFO) << msg << dt << " ms"


    // leading dimensions
    int ldA = L * a;
    int ldB = L * a;
    int ldH = G;
    int ldS = G;

    printf("[L] %d\n", L);
    printf("[A] %d\n", a);
    printf("[G] %d\n", G);
    //#ifdef GPU
    //read_clock(&t0);
    //if( myCUDAwrapper_init( tile_size ) )
    //exit(-1);
    //read_clock(&t1);
    //elapsed_us = elapsed_time( &t0, &t1 );
    //printf("[CUDA init] %.2f secs\n", elapsed_us*1e-6);
    //#endif

    //TIC;
    // H = A* T_AB B + B* T_AB* A + B* T_BB B
    read_clock(&t0);
    {
        read_clock(&t0op);
        cx_double *A_bkp; //, *B_bkp;
        if (method == "Backup") {
            // A_bkp = A
            A_bkp = (cx_double *) malloc( L * a * G * sizeof(cx_double) );
            memcpy( A_bkp, A, L * a * G * sizeof(cx_double) );

            // B_bkp = B
            //B_bkp = (cx_double *) malloc( L * a * G * sizeof(cx_double) );
            //memcpy( B_bkp, B, L * a * G * sizeof(cx_double) );
        }
        read_clock(&t1);
        elapsed_us = elapsed_time( &t0, &t1 );
        printf("[A,B backup] %.2f secs\n", elapsed_us*1e-6);

        // X stored in A
        X = A;
        int ldX = ldA;
        int X_height = 0;

#if 1
        // init Tmp
        cx_double *Tmp = (cx_double *) malloc ( L * G * sizeof(cx_double) );
        int ldTmp = L;

        printf("MKL_THREADS: %d\n", mkl_get_max_threads());

        read_clock(&t0op);
        for (size_t atom = 0; atom < a; atom++) {
            //printf("Atom %d\n", atom);
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
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[Atoms 1] %.2f secs\n", elapsed_us*1e-6);
#endif

#if 0
        read_clock(&t0op);
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

            zlacpy_("A", &L, &G, Tmp, &L, X+atom*L, &ldX);
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
#if defined GPU || defined HS_CUDA_BLAS
        gpu_zher2k_("L", "C", &G, &X_height, &c1, B, &ldB, X, &ldX, &zero, H, &ldH);
#else
        zher2k_("L", "C", &G, &X_height, &c1, B, &ldB, X, &ldX, &c0, H, &ldH);
#endif
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[H1: zher2k] %.2f secs\n", elapsed_us*1e-6);
        TOC("H zher2k", 8 * G * 1.0 * G * 1.0 * X_height, 1.0e9, "GFLOPS/s");

        // restore A and B
        read_clock(&t0op);
        if ( method == "Backup" ) {
            printf("Restoring backup\n");
            fflush(stdout);
            memcpy( A, A_bkp, L * a * G * sizeof(cx_double) );
            //memcpy( B, B_bkp, L * a * G * sizeof(cx_double) );
            free( A_bkp );
            //free( B_bkp );
        } else // not implemented
        { 
            //  construct_AB(A, B, km, ad); 
        }
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[A,B backup] %.2f secs\n", elapsed_us*1e-6);

        //free( Tmp );
        // Tmp, A_bkp, B_bkp freed implicitly
    }
    //TOC("H = A* T_AB B + B* T_AB* A + B* T_BB B     ");

    //TIC;
    // S = A* A + B* B
    //
    // S = A* A
    TIC;
    int La = L * a;
    read_clock(&t0op);
#if defined GPU || defined HS_CUDA_BLAS
    gpu_zherk_("L", "C", &G, &La, &dONE, A, &ldA, &dZERO, S, &ldS);
#else
    //zherk_("L", "C", &G, &La, &c1, A, &ldA, &c0, S, &ldS);
    zherk_("L", "C", &G, &La, &dONE, A, &ldA, &dZERO, S, &ldS);
#endif
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[S: zherk 1] %.2f secs\n", elapsed_us*1e-6);
    TOC("S zherk 1", 4.0 * G * 1.0 * G * 1.0 * La, 1.0e9, "GFLOPS/s");

    // apply u norm to B
    read_clock(&t0op);
    apply_u_norm(L, a, G, B, l_maxs, max_lmax, u_norms );
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[S: unorm] %.2f secs\n", elapsed_us*1e-6);

    // S += B* B
    TIC;
    read_clock(&t0op);
#if defined GPU || defined HS_CUDA_BLAS
    gpu_zherk_("L", "C", &G, &La, &dONE, B, &ldB, &dONE, S, &ldS);
#else
    zherk_("L", "C", &G, &La, &dONE, B, &ldB, &dONE, S, &ldS);
#endif
    read_clock(&t1op);
    elapsed_us = elapsed_time( &t0op, &t1op );
    printf("[S: zherk 2] %.2f secs\n", elapsed_us*1e-6);
    TOC("S zherk 2", 4.0 * G * 1.0 * G * 1.0 * La, 1.0e9, "GFLOPS/s");
    //read_clock(&t1op);
    //elapsed_us = elapsed_time( &t0op, &t1op );
    //printf("[S time] %.2f secs\n", elapsed_us*1e-6);
    //TOC("S = A* A + norm(B)* norm(B)                ");

    // alive: A

    //TIC;
    // H += A* T_AA A
    read_clock(&t0op);
    {
        // buffer for Chol factorization
        cx_double *Cholbuff = (cx_double *) malloc ( L * L * sizeof(cx_double) );
        int ldCholbuff = L;

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
            zlacpy_("A", &T_AA_atom_size, &T_AA_atom_size, T_AA_atom, &ldT_AA_atom, Cholbuff, &ldCholbuff);
            zpotrf_("L", &T_AA_atom_size, Cholbuff, &ldCholbuff, &info);

            cx_double *A_atom = A + atom * L;
            if (!info) { // Chol succeeded
                 int idx = a * L - X_hpd_height - T_AA_atom_size;

                // X[idx] = A[atom]
                cx_double *X_idx = X + idx;
                zlacpy_("A", &T_AA_atom_size, &G, A_atom, &ldA, X_idx, &ldA);

                // X[idx] = L* X[idx]
                ztrmm_("L", "L", "C", "N", &T_AA_atom_size, &G, &c1, Cholbuff, &ldCholbuff, X_idx, &ldA);

                X_hpd_height += T_AA_atom_size;
                nhpd++;
            } else { // Chol failed
                int idx = X_nhpd_height;

                // A[idx] = A[atom]
                cx_double *A_idx = A + idx;
                if (A_idx != A_atom)
                    zlacpy_("A", &T_AA_atom_size, &G, A_atom, &ldA, A_idx, &ldA);

                // X[idx] = T_AA[atom] A[idx]
                cx_double *X_idx = X + idx;
                zhemm_("L", "L", &T_AA_atom_size, &G, &c1, T_AA_atom, &ldT_AA_atom, A_idx, &ldA, &c0, X_idx, &ldA);

                X_nhpd_height += T_AA_atom_size;
            }
        }
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[Atoms 2] %.2f secs\n", elapsed_us*1e-6);

        //LOG(INFO) << "#hpd: " << nhpd << " / " << a;

        printf("[X_hpd_height] %d\n", X_hpd_height);
        TIC;
        read_clock(&t0op);
        // H += X[hpd]* X[hpd]
        cx_double *X_hpd = X + a * L - X_hpd_height;
#if defined GPU || defined HS_CUDA_BLAS
        gpu_zherk_("L", "C", &G, &X_hpd_height, &dONE, X_hpd, &ldX, &dONE, H, &ldH);
#else
        //zherk_("L", "C", &G, &X_hpd_height, &c1, X_hpd, &ldX, &c1, H, &ldH);
        zherk_("L", "C", &G, &X_hpd_height, &dONE, X_hpd, &ldX, &dONE, H, &ldH);
#endif
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[H2 zherk] %.2f secs\n", elapsed_us*1e-6);
        TOC("H zherk", 4 * G * 1.0 * G * X_hpd_height, 1.0e9, "GFLOPS/s");

        // H += A[nhpd]* X[nhpd]
        printf("[X_nhpd_height] %d\n", X_nhpd_height);
        TIC;
        read_clock(&t0op);
        cx_double *A_nhpd = A;
        cx_double *X_nhpd = X;
#if defined GPU || defined HS_CUDA_BLAS
        //gpu_zgemm_("C", "N", &G, &G, &X_nhpd_height, &c1, A_nhpd, &ldA, X_nhpd, &ldX, &c1, H, &ldH);
        gpu_zherkx_("L", "C", &G, &X_nhpd_height, &c1, A_nhpd, &ldA, X_nhpd, &ldX, &dONE, H, &ldH);
#else
        zgemm_("C", "N", &G, &G, &X_nhpd_height, &c1, A_nhpd, &ldA, X_nhpd, &ldX, &c1, H, &ldH);
#endif
        read_clock(&t1op);
        elapsed_us = elapsed_time( &t0op, &t1op );
        printf("[H2 zgemm] %.2f secs\n", elapsed_us*1e-6);
        TOC("H zgemm", 8 * G * 1.0 * G * X_nhpd_height, 1.0e9, "GFLOPS/s");
    }
    //TOC("H += A* T_AA A                             ");
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[Compute time] %.2f secs\n", elapsed_us*1e-6);
    
    // timings
    //tic = t_start;
    //TOC("total                                      ");
    //return dt;
    //#ifdef GPU
    //myCUDAwrapper_finalize();
    //#endif
}

