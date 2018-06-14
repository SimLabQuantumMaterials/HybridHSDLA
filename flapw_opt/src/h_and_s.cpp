/** \file 
 * Driver program to read dumped data from Fleur and generate the overlap and
 * Hamiltonian matrix.
 * TODO: Muffin-Tin part works, interstitial missing
 */

//#define ARMA_32BIT_WORD


#include "configured_input.h"
#include "core.h"
#include "util.h"
#include "types.h"
#include "FleurData.h"
#include "HDF5Store.h"
#include "mkl_service.h"
#include <iostream>
#include <armadillo>
#include <complex>
#include <cassert>
#include <chrono>
#include <cmath>

#include "logging.h"

extern "C" {
#include "timing.h"
}
#include "h_and_s_stripped.h"
#include "types4stripped.h"
#include "blas.h"

#ifdef GPU
    #include <cuda_runtime.h>
#endif
#if defined GPU || defined HS_HSTREAMS || defined HS_CUDA
    #include "wrappers.h"
#endif

using namespace flapw;

// sample data used if none given on command line
const extern char* DEFAULT_DUMP_DIR;

/*
HDF5Store *pdump = NULL;
#ifdef NDEBUG
#define DUMP(x)
#else
#define DUMP(x) {pdump->add_arma((x), #x); pdump->flush();}
#endif
*/

/** \brief Generate the overlap matrix
 *
 * "Naive" idea:
 * Generate full multi-index AB matrix
 * -> calculate S from A(k), B(k), as
 *  S_ij(k) = \sum_{m, L} A(k)*_{L, m, i} A(k)_{L, m, j}
 *                      + B(k)*_{L, m, i} B(k)_{L, m, j} * dE_u_norm(L, m)
 *
 * Implemented as a HERK-Blas3-call, which should give much higher performance.
 */
arma::cx_mat contract_AB_overlap(
        const arma::cx_cube& A,
        const arma::cx_cube& B,
        const match_field& mf)
{
    const arma::blas_int n = A.n_slices;
    // "reshape" tensor, treat first two dimensions as one
    const arma::blas_int size_rc = A.n_rows * A.n_cols;
    arma::cx_mat S(n, n);
    // first sum, herk(A) = A^h * A
    double alpha = 1.;
    double beta = 0.;
    arma::blas::herk("L", "C", &n, &size_rc, &alpha, A.memptr(), &size_rc, &beta, S.memptr(), &n);
    assert( (B.n_slices == n) && (B.n_rows*B.n_cols == size_rc) );
    // second sum, herk(B)*dE_u_norm, sqrt(dE_u_norm) implicitly contained in B
    beta = 1.; // add to first sum
    arma::blas::herk("L", "C", &n, &size_rc, &alpha, B.memptr(), &size_rc, &beta, S.memptr(), &n);
    return S;
}

/** \brief naive (non-BLAS) version, directly implements the formula
 *
 * \see contract_AB_overlap for details
 */
arma::cx_mat contract_AB_overlap_naive(
            const arma::cx_cube& A,
            const arma::cx_cube& B,
            const match_field& mf)
{
    int n = A.n_slices;
    int max_lmax = mf.n_rows - 1;
    arma::cx_mat S(n, n);
    for (size_t G = 0; G < n; ++G) {
        for (size_t G_prime = 0; G_prime < n; ++G_prime) {
            S(G_prime, G) = 0;
            for (size_t at = 0; at < A.n_cols; ++at) {
                for (int l = 0; l <= max_lmax; ++l) {
                    auto dE_u_norm = mf(l, at).dE_u_norm;
                    for (int m = -l; m <= l; ++m) {
                        int lm = l*(l+1) + m;
                        S(G_prime, G) +=
                            std::conj(A(lm, at, G_prime)) * A(lm, at, G)
                          + std::conj(B(lm, at, G_prime))
                                    * B(lm, at, G) * dE_u_norm;
                    }
                }
            }
        }
    }
    return S;
}
/** \brief Construct the Hamilton matrix from AB tensors and T - direct sum
 *
 * Given the matching coefficients and the four T-matrices, calculate the
 * complete Hamiltonian by executing the given sum element-by-element.
 *
 * \param       A   3-tensor of function matching coefficients, size (L,a,n)
 * \param       B   similar to A, but matching coefficients of the derivative
 * \param       T   collects the 4 T-matrices of size (L,L)
 * \returns the complete Hamiltonian H of size n^2
 */
arma::cx_mat contract_AB_hamilton_naive(
                const arma::cx_cube& A,
                const arma::cx_cube& B,
                const AtomData::TMatVec& T)
{
    assert(T.size() == A.n_cols); // one entry per atom
    int n = A.n_slices;
    arma::cx_mat H(n, n);
    for (size_t G = 0; G < n; ++G) {
      for (size_t G_prime = 0; G_prime < n; ++G_prime) {
        H(G, G_prime) = 0;
        for (size_t at = 0; at < A.n_cols; ++at) {
          for (size_t lm = 0; lm < T[at].n; ++lm) {
            for (size_t lm_prime = 0; lm_prime < T[at].n; ++lm_prime) {
                H(G, G_prime) +=  std::conj(A(lm_prime, at, G))
                                * T[at].aa(lm_prime, lm)
                                * A(lm, at, G_prime);
                H(G, G_prime) +=  std::conj(B(lm_prime, at, G))
                                * T[at].bb(lm_prime, lm)
                                * B(lm, at, G_prime);
                H(G, G_prime) +=  std::conj(A(lm_prime, at, G))
                                * T[at].ab(lm_prime, lm)
                                * B(lm, at, G_prime);
                H(G, G_prime) +=  std::conj(B(lm_prime, at, G))
                                * T[at].ba(lm_prime, lm)
                                * A(lm, at, G_prime);
            }
          }
        }
      }
    }
#ifndef NDEBUG
    if (!is_hermitian(H)) {
            LOG(ERROR) << "H is not Hermitian";
    }
#endif
    return H;
}

/** \brief Execute a single part of the Hamilton summation
 *
 * Calculate a contraction of the form
 * \f[
 * H_{lm} = \sum_{ik} U^*_{il} T_{ik} V_{km}
 * \f]
 */
void hamilton_kernel(const arma::cx_cube& U,
                     const arma::cx_cube& V,
                     const arma::cx_mat& T,
                     const int a_offset,
                     const int L_nonsph,
                     const int L_diag,
                     arma::cx_mat& H)
{
    int Lmax = U.n_rows, num_a = U.n_cols, num_g = U.n_slices;
    cx_double c_one = 1., c_zero = 0.;
    // actual GEMM may be smaller than Lmax:
    int num_L = T.n_rows;
    int LDUV = Lmax*num_a;
    int el_offset = Lmax*a_offset;
    arma::cx_mat tmp(num_L, num_g); // storage for first GEMM
    ///////////////////////////////
    // first GEMM with the middle and last term, tmp = T * V
    ///////////////////////////////
    LOG(DEBUG) << "GEMM tmp=T*V: m=" << num_L << ", n=" << num_g << ", k=" << num_L
               << "\nLDA=" << num_L << ", el_offset=" << el_offset << "(V.memptr()=" << V.memptr()
               << "; V_ptr +el_off=" << V.memptr()+el_offset << "), LDB=" << LDUV;
    // dense, nonspherical part
    arma::blas::gemm("N", "N", &L_nonsph, &num_g, &L_nonsph,
                     &c_one, T.memptr(), &num_L,
                     V.memptr() + el_offset, &LDUV,
                     &c_zero, tmp.memptr(), &num_L);
    // diagonal, spherical-only part
    if (L_diag > 0) {
        LOG(DEBUG) << "Diagonal extra: size = " << L_diag;
        // extract spherical-only diagonal
        arma::cx_vec T_diag = T.diag();
        T_diag = T_diag.rows(L_nonsph, num_L-1);
        for (int g = 0; g < num_g; ++g) {
            // for each column 'g', set lower rows
            tmp(L_nonsph, g, arma::size(L_diag, 1)) = T_diag % V.slice(g)(L_nonsph, a_offset, arma::size(L_diag, 1));
            //(L_nonsph, a_offset, g, arma::size(L_diag, 1, 1));
            //tmp(L_nonsph, g, arma::size(L_diag, 1)) %= T_diag;
        }
    }

    LOG(DEBUG) << "GEMM H=U*tmp: m=" << num_g << ", n=" << num_g << ", k=" << num_L
               << "\nLDA=" << LDUV << ", el_offset=" << el_offset << "(U.memptr()=" << U.memptr()
               << "; U_ptr +el_off=" << U.memptr()+el_offset << "), LDB=" << num_L;
    // second GEMM: remaining factor with result, U^H * tmp
    arma::blas::gemm("C", "N", &num_g, &num_g, &num_L,
                     &c_one, U.memptr() + el_offset, &LDUV,
                     tmp.memptr(), &num_L,
                     &c_one, H.memptr(), &num_g);
}


/** \brief Construct the Hamilton matrix from AB tensors and T - BLAS
 *
 * Using BLAS, contract the A,B,T matrices with multiple GEMM calls. For parameter
 * description and general information, see contract_AB_hamilton_naive().
 *
 * TODO: huge workspace matrix needed, blocking approach better?
 */
arma::cx_mat contract_AB_hamilton(
                const arma::cx_cube& A,
                const arma::cx_cube& B,
                const AtomData::TMatVec& T)
{
    assert(T.size() == A.n_cols); // one entry per atom
    arma::cx_mat H(A.n_slices, A.n_slices, arma::fill::zeros);
    LOG(DEBUG) << "Contracting H for " << A.n_cols << " atoms.";
    // handle AB term first, so we can do without explicit multiplication for BA
    for (size_t at = 0; at < A.n_cols; ++at) {
        hamilton_kernel(A, B, T[at].ab, at, T[at].L_nonspherical, T[at].L_diag, H);
    }
    // works: A^H T_AB B == [B^H T_BA A]^H, since T_AB = T_BA^H by construction
    H += H.t(); // handle BA term
    for (size_t at = 0; at < A.n_cols; ++at) {
        LOG(DEBUG) << "... atom " << at;
        hamilton_kernel(A, A, T[at].aa, at, T[at].L_nonspherical, T[at].L_diag, H);
        hamilton_kernel(B, B, T[at].bb, at, T[at].L_nonspherical, T[at].L_diag, H);
#ifndef NDEBUG
        if (!is_hermitian(H)) {
            LOG(ERROR) << "H is not Hermitian: atom " << at;
        }
#endif
    }
    return H;
}

std::string t_pretty(double dt) 
{
    std::ostringstream os;
    int h,m, s;
    s = dt/1000; //secs
    h = s / 3600;
    m = (s % 3600)/60;
    dt = std::fmod(dt/1000, 60);
    os << "(" << h << ":" << m << ":" << dt << ")";
    return os.str();
}

enum class GenMethod { compute, store };

double create_HS(const arma::vec3& kpt,
                 const double kmax,
                 const arma::mat33& bmat,
                 const arma::imat& g_mesh,
                 const AtomData& ad,
                 const GenMethod method,
                 arma::cx_mat& H_out,
                 arma::cx_mat& S_out)
{
    // timing stuff
    using Clock = std::chrono::system_clock;
    auto diff_as_msecs = [](Clock::time_point t) {
        auto d = Clock::now() - t;
        return std::chrono::duration_cast<
            std::chrono::milliseconds>(d).count();
    };
    auto t_start = Clock::now();
    decltype(t_start) tic;
    double dt;
#define TIC tic = Clock::now()
#define TOC(msg) dt = diff_as_msecs(tic); LOG(INFO) << msg << dt << " ms"

    auto km = create_kmesh(kpt, kmax, bmat, g_mesh);

    TIC;
    arma::cx_cube A, B;
    construct_AB(A, B, km, ad);
    TOC("A, B = construct()                         ");

    auto u_match = ad.get_u_match();
    const auto& tmats = ad.get_tmats_handle();

    TIC;
    H_out = contract_AB_hamilton(A, B, tmats);
    TOC("H = ...                                    ");

    TIC;
    ad.apply_u_norm(B);
    S_out = contract_AB_overlap(A, B, u_match);
    TOC("S = ...                                    ");

    // timings
    tic = t_start;
    TOC("total                                      ");
    return dt;
}

// declare some BLAS/LAPACK signatures
#if defined WITH_BLAS
extern "C" {
void zhemm(const char *, const char *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);
void zher2k(const char *, const char *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);
void zlacpy(const char *, const int *, const int *, const cx_double *, const int *, cx_double *, const int *);
void zpotrf(const char *, const int *, cx_double *, const int *, int *);
void ztrmm(const char *, const char *, const char *, const char *, const int *, const int *, const cx_double *, const cx_double *, const int *, cx_double *, const int *);
void zherk(const char *, const char *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);
void zgemm(const char *, const char *, const int *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);
}

/*
 *                                  alive   needs   destroys
 * A, B = construct()
 * H = A* T B + B* T A + B* T B     A, B    A, B    A
 * A = (re)construct()              B
 * S = A* A + B* B                  A, B    A, B    B
 * H += A* TA                       A,      A       A, X
 */
double create_HS2(const arma::vec3& kpt,
                  const double kmax,
                  const arma::mat33& bmat,
                  const arma::imat& g_mesh,
                  const AtomData& ad,
                  const GenMethod method,
                  arma::cx_mat& H,
                  arma::cx_mat& S)
{
    LOG(INFO) << "CPU BLAS Version - ELMAR";
    // constants
    const cx_double c0(0, 0);
    const cx_double c1(1, 0);
    const cx_double cp5(.5, 0);

    // timing stuff
    using Clock = std::chrono::system_clock;
    auto diff_as_msecs = [](Clock::time_point t) {
        auto d = Clock::now() - t;
        return std::chrono::duration_cast<
            std::chrono::milliseconds>(d).count();
    };
    auto t_start = Clock::now();
    decltype(t_start) tic;
    double dt;
#define TIC tic = Clock::now()
#define TOC(msg) dt = diff_as_msecs(tic); LOG(INFO) << msg << dt << " ms"

    decltype(t_start) tic2;
    double dt2;
    double gflops;
#define TIC_NESTED tic2 = Clock::now()
#define TOC_NESTED(msg) dt2 = diff_as_msecs(tic2); LOG(INFO) << msg << dt2 << " ms" << "    " << gflops/(dt2*1e-3)

    // get some pysics stuff
    auto km = create_kmesh(kpt, kmax, bmat, g_mesh);

    // get Ts
    const auto &T = ad.get_tmats_handle();

    TIC;
    // create A, B
    arma::cx_cube A, B;
    construct_AB(A, B, km, ad);
    TOC("A, B = construct()                         ");

    // extract sizes
    const int L = A.n_rows;
    const int a = A.n_cols;
    const int G = A.n_slices;

    LOG(INFO) << "L = " << L << "; a = " << a << "; G = " << G;

    // leading dimensions
    const int ldA = L * a;
    const int ldB = L * a;
    const int ldH = G;
    const int ldS = G;

    // init H
    H.set_size(G, G);

    TIC;
    // H = A* T_AB B + B* T_AB* A + B* T_BB B
    {
        arma::cx_cube A_bkp, B_bkp;
        if (method == GenMethod::store) {
            // A_bkp = A
            A_bkp = arma::cx_cube(A);

            // B_bkp = B
            B_bkp = arma::cx_cube(B);
        }

        // X stored in A
        arma::cx_cube &X = A;
        const int ldX = ldA;

        int X_height = 0;

        // init Tmp
        arma::cx_mat Tmp_mat(L, G);
        cx_double *Tmp = Tmp_mat.memptr();
        const int ldTmp = L;

        //mkl_set_num_threads(24);
        TIC_NESTED;
        for (size_t atom = 0; atom < a; atom++) {
            // Tmp = 0
            Tmp_mat.zeros();

            // Tmp = T_BB[atom] B[atom]
            const arma::cx_mat &T_BB_atom_mat = T[atom].bb;
            const int T_BB_atom_size = T_BB_atom_mat.n_rows;
            const int ldT_BB_atom = T_BB_atom_size;
            const cx_double *T_BB_atom = T_BB_atom_mat.memptr();
            const cx_double *B_atom = B.memptr() + atom * L;
            zhemm("L", "L", &T_BB_atom_size, &G, &cp5, T_BB_atom, &ldT_BB_atom, B_atom, &ldB, &c0, Tmp, &ldTmp);

            // Tmp += T_AB[atom]* A[atom]
            const arma::cx_mat &T_AB_atom_mat = T[atom].ab;
            const int T_AB_atom_size = T_AB_atom_mat.n_rows;
            const int ldT_AB_atom = T_AB_atom_size;
            const cx_double *T_AB_atom = T_AB_atom_mat.memptr();
            const cx_double *A_atom = A.memptr() + atom * L;
            zgemm("C", "N", &T_AB_atom_size, &G, &T_AB_atom_size, &c1, T_AB_atom, &ldT_AB_atom, A_atom, &ldA, &c1, Tmp, &ldTmp);

            // X[X_height] = Tmp
            const int Tmp_size = (T_AB_atom_size > T_BB_atom_size) ? T_AB_atom_size : T_BB_atom_size;
            cx_double *X_X_height = X.memptr() + X_height;
            zlacpy("A", &Tmp_size, &G, Tmp, &ldTmp, X_X_height, &ldX);

            if (X_height != atom * L) {
                // Tmp = B_atom
                const cx_double *B_atom = B.memptr() + atom * L;
                zlacpy("A", &Tmp_size, &G, B_atom, &ldB, Tmp, &ldTmp);

                // B[X_height] = Tmp
                cx_double *B_X_height = B.memptr() + X_height;
                zlacpy("A", &Tmp_size, &G, Tmp, &ldTmp, B_X_height, &ldTmp);
            }

            X_height += Tmp_size;
        }
        gflops = 0.0;
        TOC_NESTED("Loop 1: ");

        // H = B* X + X* B
        const cx_double *B_all = B.memptr();
        const cx_double *X_all = X.memptr();
        cx_double *H_all = H.memptr();

        printf("X height: %d\n", X_height);
        printf("G: %d\n", G);
        gflops = (double)G * X_height * G * 8 * 1e-9;
        printf("GFs: %f\n", gflops);
        TIC_NESTED;
        //printf("MKL %d\n", mkl_get_max_threads());
        zher2k("L", "C", &G, &X_height, &c1, B_all, &ldB, X_all, &ldX, &c0, H_all, &ldH);
        TOC_NESTED("ZHERK: ");

        // restore A and B
        if (method == GenMethod::store) {
            A = arma::cx_cube(A_bkp);
            B = arma::cx_cube(B_bkp);
        } else
            construct_AB(A, B, km, ad); 

        // Tmp, A_bkp, B_bkp freed implicitly
    }
    TOC("H = A* T_AB B + B* T_AB* A + B* T_BB B     ");

    // init S
    S.set_size(G, G);

    TIC;
    // S = A* A + B* B
    {
        // S = A* A
        const int La = L * a;
        const cx_double *A_all = A.memptr();
        cx_double *S_all = S.memptr();
        //gflops = G * La * G * 4 * 1e-9;
        //TIC_NESTED;
        zherk("L", "C", &G, &La, &c1, A_all, &ldA, &c0, S_all, &ldS);
        //TOC_NESTED("ZHERK: ");

        // apply u norm to B
        ad.apply_u_norm(B);

        // S += B* B
        const cx_double *B_all = B.memptr();
        zherk("L", "C", &G, &La, &c1, B_all, &ldB, &c1, S_all, &ldS);
    }
    TOC("S = A* A + norm(B)* norm(B)                ");

    // alive: A

    TIC;
    // H += A* T_AA A
    {
        // buffer for Chol factorization
        arma::cx_mat Cholbuff_mat(L, L);
        cx_double *Cholbuff = Cholbuff_mat.memptr();
        const int ldCholbuff = L;

        // X stored in B
        arma::cx_cube &X = B;
        const int ldX = ldB;

        int nhpd = 0;
        int X_nhpd_height = 0;
        int X_hpd_height = 0;

        for (size_t atom = 0; atom < a; atom++) {
            // Cholbuff = Chol(T_AA[atom])
            const arma::cx_mat &T_AA_atom_mat = T[atom].aa;
            const int T_AA_atom_size = T_AA_atom_mat.n_rows;
            const int ldT_AA_atom = T_AA_atom_size;
            const cx_double *T_AA_atom = T_AA_atom_mat.memptr();
            int info;
            zlacpy("A", &T_AA_atom_size, &T_AA_atom_size, T_AA_atom, &ldT_AA_atom, Cholbuff, &ldCholbuff);
            zpotrf("L", &T_AA_atom_size, Cholbuff, &ldCholbuff, &info);

            cx_double *A_atom = A.memptr() + atom * L;
            if (!info) { // Chol succeeded
                const int idx = a * L - X_hpd_height - T_AA_atom_size;

                // X[idx] = A[atom]
                cx_double *X_idx = X.memptr() + idx;
                zlacpy("A", &T_AA_atom_size, &G, A_atom, &ldA, X_idx, &ldA);

                // X[idx] = L* X[idx]
                ztrmm("L", "L", "C", "N", &T_AA_atom_size, &G, &c1, Cholbuff, &ldCholbuff, X_idx, &ldA);

                X_hpd_height += T_AA_atom_size;
                nhpd ++;
            } else { // Chol failed
                const int idx = X_nhpd_height;

                // A[idx] = A[atom]
                cx_double *A_idx = A.memptr() + idx;
                if (A_idx != A_atom)
                    zlacpy("A", &T_AA_atom_size, &G, A_atom, &ldA, A_idx, &ldA);

                // X[idx] = T_AA[atom] A[idx]
                cx_double *X_idx = X.memptr() + idx;
                zhemm("L", "L", &T_AA_atom_size, &G, &c1, T_AA_atom, &ldT_AA_atom, A_idx, &ldA, &c0, X_idx, &ldA);

                X_nhpd_height += T_AA_atom_size;
            }
        }

        LOG(INFO) << "#hpd: " << nhpd << " / " << a;

        // H += X[hpd]* X[hpd]
        const cx_double *X_hpd = X.memptr() + a * L - X_hpd_height;
        cx_double *H_all = H.memptr();
        zherk("L", "C", &G, &X_hpd_height, &c1, X_hpd, &ldX, &c1, H_all, &ldH);

        // H += A[nhpd]* X[nhpd]
        const cx_double *A_nhpd = A.memptr();
        const cx_double *X_nhpd = X.memptr();
        zgemm("C", "N", &G, &G, &X_nhpd_height, &c1, A_nhpd, &ldA, X_nhpd, &ldX, &c1, H_all, &ldH);
    }
    TOC("H += A* T_AA A                             ");
    
    // timings
    tic = t_start;
    TOC("total                                      ");
    return dt;
}

#elif defined STRIPPED
double create_HS2(const arma::vec3& kpt,
                  const double kmax,
                  const arma::mat33& bmat,
                  const arma::imat& g_mesh,
                  const AtomData& ad,
                  const GenMethod method,
                  arma::cx_mat& H,
                  arma::cx_mat& S)
{
    LOG(INFO) << "Refactored Version - HybridTCon";
    // constants
    const cx_double c0(0, 0);
    const cx_double c1(1, 0);
    const cx_double cp5(.5, 0);

    // timing stuff
    using Clock = std::chrono::system_clock;
    auto diff_as_msecs = [](Clock::time_point t) {
        auto d = Clock::now() - t;
        return std::chrono::duration_cast<
            std::chrono::milliseconds>(d).count();
    };
    auto t_start = Clock::now();
    decltype(t_start) tic;
    double dt;
#define TIC tic = Clock::now()
#define TOC(msg) dt = diff_as_msecs(tic); LOG(INFO) << msg << dt << " ms"

    // more timings
    struct timeval t0, t1;
    long int elapsed_us;

    printf("[START]\n");

    // get some pysics stuff
    read_clock(&t0);
    auto km = create_kmesh(kpt, kmax, bmat, g_mesh);
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[Create kmesh] %.3f secs\n", elapsed_us*1e-6);

    // get Ts
    const auto &T = ad.get_tmats_handle();

    // [DIEGO] extract sizes
    const int L = (ad.max_lmax+1) * (ad.max_lmax+1);
    const int a = ad.num;
    const int G = km.num;
    LOG(INFO) << "L = " << L << "; a = " << a << "; G = " << G;

    // [DIEGO] Reduce memory footprint of tests by pinning A/B directly
    read_clock(&t0);
    cx_double *A_pinned;
    cx_double *B_pinned;
    cx_double *H_pinned;
    cx_double *S_pinned;

    // Allocate pinned workspace for computing A and B
    cx_double *W_pinned;
#ifdef GPU
    if( cudaHostAlloc( (void **)&A_pinned, (size_t)L * a * G * sizeof(cx_double), 0 ) != cudaSuccess )
        printf("[ERROR] Couldn't allocate pinned memory for A\n");
    if( cudaHostAlloc( (void **)&B_pinned, (size_t)L * a * G * sizeof(cx_double), 0 ) != cudaSuccess )
        printf("[ERROR] Couldn't allocate pinned memory for B\n");
    if( cudaHostAlloc( (void **)&W_pinned, (size_t)L * a * G * sizeof(cx_double), 0 ) != cudaSuccess )
        printf("[ERROR] Couldn't allocate pinned memory for W\n");
#else
    A_pinned = (cx_double *) malloc ((size_t)L * G * a * sizeof(cx_double));
    B_pinned = (cx_double *) malloc ((size_t)L * G * a * sizeof(cx_double));
    W_pinned = (cx_double *) malloc ((size_t)L * G * a * sizeof(cx_double));
#endif
    arma::cx_cube A(A_pinned, L, a, G, false, true),
                  B(B_pinned, L, a, G, false, true);
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[CUDA HostAlloc 2] %.2f secs\n", elapsed_us*1e-6);


    TIC;
    // create A, B
    //arma::cx_cube A, B;
    //
    read_clock(&t0);
    construct_AB(A, B, km, ad);
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[Construct A,B] %.3f secs\n", elapsed_us*1e-6);
    //
    TOC("A, B = construct()                         ");


    // [DIEGO]
    // Prepare all data for a regular non-armadillo function
    //
    read_clock(&t0);
    // Buffers
    // A, B
    cx_double *A_ptr = A.memptr();
    cx_double *B_ptr = B.memptr();
    // H, S
    H.set_size(G, G);
    cx_double *H_ptr = H.memptr();
    S.set_size(G, G);
    cx_double *S_ptr = S.memptr();
    // A,B creation/storage method
    std::string method_str = "Backup";
    //
    assert( a == T.num );
    // l_max per atom (first entry of the tuple)
    int *lmaxs = (int *) malloc( a * sizeof(int) );
    for ( int at = 0; at < a; at++ )
        lmaxs[ at ] = ad.get_lmax(at).first;
    int max_lmax = ad.max_lmax;
    // Collect sqrt of u_norms
    double *sqrt_u_norms = (double *) std::calloc( (max_lmax + 1) * a, sizeof(double) );
    for ( int at = 0; at < a; at++ )
    {
        int lmax = lmaxs[ at ];
        for ( int l = 0; l <= lmax; l++ )
            sqrt_u_norms[ at*(max_lmax+1) + l ] = sqrt(ad.u_match(l, at).dE_u_norm);
    }
    if ( lmaxs == NULL || sqrt_u_norms == NULL )
    {
        printf("Cannot allocate lmaxs || sqrt_u_norms\n");
        exit(-2);
    }
    // TMats
    myTMat *myTMats = (myTMat *) malloc( a * sizeof(myTMat) );
    if ( myTMats == NULL )
    {
        printf("Cannot allocate myTMat\n");
        exit(-2);
    }
    for ( int atom = 0; atom < a; atom++ )
    {
        // AA
        const arma::cx_mat &T_AA_atom_mat = T[atom].aa;
        const int T_AA_atom_size = T_AA_atom_mat.n_rows;
        const cx_double *T_AA_atom = T_AA_atom_mat.memptr();
        const arma::cx_mat &T_AB_atom_mat = T[atom].ab;
        const int T_AB_atom_size = T_AB_atom_mat.n_rows;
        const cx_double *T_AB_atom = T_AB_atom_mat.memptr();
        const arma::cx_mat &T_BB_atom_mat = T[atom].bb;
        const int T_BB_atom_size = T_BB_atom_mat.n_rows;
        const cx_double *T_BB_atom = T_BB_atom_mat.memptr();
        myTMats[ atom ].set( T_AA_atom_size, T_AB_atom_size, T_BB_atom_size,
                             T_AA_atom, T_AB_atom, T_BB_atom );
    }
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[Preparing data] %.3f secs\n", elapsed_us*1e-6);
    // [/ DIEGO]

#ifdef DUMP_FLAPW_DATA
    // [DUMP]
    char *dump_dir_str, *end;
    char dump_file_path[1024];
    dump_dir_str = getenv( "FLAPW_DUMP_DIR" );
    
    if ( dump_dir_str == NULL )
    {
        fprintf("[ERROR] Environment variable FLAPW_DUMP_DIR not set! Required for data dumping.");
        exit(-1);
    }

    snprintf( dump_file_path, 1024, "%s/data.in", dump_dir_str );
    FILE *fp = fopen( dump_file_path, "wb" );
    if ( fp == NULL )
    {
        printf("Cannot open %s\n", dump_file_path);
        exit(-3);
    }
    // Write problem dimensions: L, a, G
    fwrite( &L, sizeof(int), 1, fp );
    fwrite( &a, sizeof(int), 1, fp );
    fwrite( &G, sizeof(int), 1, fp );
    // Write A and B buffers
    fwrite( A_ptr, sizeof(cx_double), L*a*G, fp );
    fwrite( B_ptr, sizeof(cx_double), L*a*G, fp );
    // Write T_xx matrices
    for ( int atom = 0; atom < a; atom++ )
    {
        int aa_size = myTMats[ atom ].aa_size;
        fwrite( &aa_size, sizeof(int), 1, fp );
        fwrite( myTMats[ atom ].aa, sizeof(cx_double), aa_size*aa_size, fp );
        int ab_size = myTMats[ atom ].ab_size;
        fwrite( &ab_size, sizeof(int), 1, fp );
        fwrite( myTMats[ atom ].ab, sizeof(cx_double), ab_size*ab_size, fp );
        int bb_size = myTMats[ atom ].bb_size;
        fwrite( &bb_size, sizeof(int), 1, fp );
        fwrite( myTMats[ atom ].bb, sizeof(cx_double), bb_size*bb_size, fp );
    }
    // Write buffers to compute norms
    fwrite( lmaxs, sizeof(int), a, fp );
    fwrite( &max_lmax, sizeof(int), 1, fp );
    fwrite( sqrt_u_norms, sizeof(double), (max_lmax+1)*a, fp );
    // Done
    fclose( fp );
    // [/ DUMP]
#endif

    read_clock(&t0);
#if defined GPU || defined HS_HSTREAMS || defined HS_CUDA
    if( dev_init( ) )
        exit(-1);
#endif
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[CUDA Init] %.3f secs\n", elapsed_us*1e-6);

    read_clock(&t0);
#ifdef GPU
    printf("Allocating pinned memory\n");
    //if( cudaHostAlloc( (void **)&A_pinned, (size_t)L * a * G * sizeof(cx_double), 0 ) != cudaSuccess )
    //printf("[ERROR] Couldn't allocate pinned memory for A\n");
    //if( cudaHostAlloc( (void **)&B_pinned, (size_t)L * a * G * sizeof(cx_double), 0 ) != cudaSuccess )
    //printf("[ERROR] Couldn't allocate pinned memory for B\n");
    if( cudaHostAlloc( (void **)&H_pinned, (size_t)G * G * sizeof(cx_double), 0 ) != cudaSuccess )
        printf("[ERROR] Couldn't allocate pinned memory for H\n");
    if( cudaHostAlloc( (void **)&S_pinned, (size_t)G * G * sizeof(cx_double), 0 ) != cudaSuccess )
        printf("[ERROR] Couldn't allocate pinned memory for S\n");

    //memcpy( A_pinned, A_ptr, (size_t)L * a * G * sizeof(cx_double) );
    //memcpy( B_pinned, B_ptr, (size_t)L * a * G * sizeof(cx_double) );
    //memcpy( H_pinned, H_ptr, (size_t)G * G * sizeof(cx_double) );
    //memcpy( S_pinned, S_ptr, (size_t)G * G * sizeof(cx_double) );
#else
    A_pinned = A_ptr;
    B_pinned = B_ptr;
    //A_pinned = (cx_double *) malloc ((size_t)L*G*a*sizeof(cx_double));
    //B_pinned = (cx_double *) malloc ((size_t)L*G*a*sizeof(cx_double));
    H_pinned = H_ptr;
    S_pinned = S_ptr;
#endif
    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[CUDA HostAlloc] %.3f secs\n", elapsed_us*1e-6);

    //create_HS2_stripped_v2( L, a, G, A_ptr, A_pinned, B_ptr, B_pinned, 
    create_HS2_stripped_v2( L, a, G, A_pinned, B_pinned, H_pinned, S_pinned,
                         myTMats, lmaxs, max_lmax, sqrt_u_norms,
                         method_str, W_pinned );

    read_clock(&t0);
#ifdef GPU
    //read_clock(&t0);
    //memcpy( A_ptr, A_pinned, (size_t)L * a * G * sizeof(cx_double) );
    //memcpy( B_ptr, B_pinned, (size_t)L * a * G * sizeof(cx_double) );
    memcpy( H_ptr, H_pinned, (size_t)G * G * sizeof(cx_double) );
    memcpy( S_ptr, S_pinned, (size_t)G * G * sizeof(cx_double) );

    cudaFreeHost( H_pinned );
    cudaFreeHost( S_pinned );
    cudaFreeHost( A_pinned );
    cudaFreeHost( B_pinned );
    cudaFreeHost( W_pinned );
    //read_clock(&t1);
    //elapsed_us = elapsed_time( &t0, &t1 );
    //printf("[Pinned to local] %.3f secs\n", elapsed_us*1e-6);

    //dev_finalize();
#else
    //free(A_pinned);
    //free(B_pinned);
    free(W_pinned);
#endif
    //dev_finalize();

    read_clock(&t1);
    elapsed_us = elapsed_time( &t0, &t1 );
    printf("[CUDA Free & Finalize] %.3f secs\n", elapsed_us*1e-6);

#ifdef DUMP_FLAPW_DATA
    // [DUMP Results]
    snprintf( dump_file_path, 1024, "%s/data.out", dump_dir_str );
    FILE *fp = fopen( dump_file_path, "wb" );
    if ( fp == NULL )
    {
        printf("Cannot open %s\n", dump_file_path);
        exit(-3);
    }
    fwrite( &G, sizeof(int), 1, fp );
    fwrite( H_ptr, sizeof(cx_double), G*G, fp );
    fwrite( S_ptr, sizeof(cx_double), G*G, fp );
    // Done
    fclose( fp );
    // [/ DUMP Results]
#endif

    free( lmaxs );
    free( sqrt_u_norms );
    free( myTMats );

    // timings
    tic = t_start;
    TOC("total                                      ");

    printf("[TOTAL] %.3f secs\n", dt*1e-3);

    printf("[END]\n");

    return dt;
}
#endif

/* driver function */
void generate_matrices(const std::string& dump_dir,
                      int k_ind_max,
                      GenMethod method)
{
    const std::string hdf_out_path = "flapw_out.h5";
    FleurData fdata(dump_dir);
    HDF5Store hdf_out(hdf_out_path);
//     pdump = &hdf_out;
    k_ind_max = std::min(k_ind_max, fdata.get_num_kpts());

    double kmax = 5.;
    auto bmat = fdata.get_reciprocal_basis();
    auto ad = fdata.populate_atomdata();

    ad.generate_matchdata();
    ad.augment_tmats();
    arma::vec times(k_ind_max);
    for (int k_ind = 0; k_ind < k_ind_max; ++k_ind) {
        LOG(INFO) << "Starting k-point " << k_ind;
        auto g_mesh = fdata.get_g_mesh(k_ind);
        auto k_pt = fdata.get_kpoint(k_ind);
        arma::cx_mat H_out, S_out;
        times(k_ind) = create_HS2(k_pt, kmax, bmat, g_mesh, ad, method, H_out, S_out);

        continue;

        hdf_out.add_hsmat(H_out, S_out, k_ind);
        arma::cx_mat H_out2, S_out2;
        create_HS(k_pt, kmax, bmat, g_mesh, ad, method, H_out2, S_out2);
        // continue;
        const auto diffH = arma::norm(arma::trimatl(H_out2 - H_out), "inf");
        const auto diffS = arma::norm(S_out2 - S_out, "inf");
        LOG(INFO) << "norm(difference in H): " << diffH;
        LOG(INFO) << "norm(difference in S): " << diffS;
    }
    // LOG(INFO) << "Data path was " << dump_dir;
    // LOG(INFO) << "H, S saved to " << hdf_out_path;
    times.save("time.info",arma::arma_ascii);
}


int main(int argc, char** argv) {
    std::string dump_path = "";
    int k_ind_max = 60;
    GenMethod m = GenMethod::store;
    if (argc > 1)
        dump_path = argv[1];
    if (argc > 2)
        k_ind_max = atoi(argv[2]);
    if (argc > 3)
        m = std::string(argv[3]) == "compute" ? GenMethod::compute : m;
    if (dump_path == "") dump_path = DEFAULT_DUMP_DIR;
    // LOG(INFO) << "Starting. Data path = " << dump_path 
    //           << ", max. k-points = " << k_ind_max;
    // int num_threads = mkl_get_max_threads();
    // LOG(INFO) << "MKL: Using " << num_threads << " threads.";
    generate_matrices(dump_path, k_ind_max, m);
}
