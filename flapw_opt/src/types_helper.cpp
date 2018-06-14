/** \file
 * Additional function to handle/create primitive simulation types.
 */
#include "types.h"
#include "binmat.h"
#include "util.h"
#include <cassert>

namespace flapw {

/** \brief Create the k-mesh for given k-point, kmax and basis vectors (bmat)
 *
 * Currently reads the dump of a k-mesh from Fleur, and supplies the kpt instead of
 * taking it as an input.
 * \param         kpt          k-point around which the mesh is constructed
 * \param[in]     kmax         radius of the k-sphere which contains all g-points
 * \param[in]     bravais_mat  matrix of basis vectors (cols) of the reciprocal space
 * \returns the k-mesh corresponding to parameter `kpt`, with some metadata, \see struct k_mesh.
 */
k_mesh create_kmesh(const arma::vec3& kpt, const double kmax,
                    const arma::mat33& bravais_mat, arma::imat g_dump)
{
    k_mesh km;
    km.g = g_dump; // TODO: generate, don't use dump
    assert(km.g.n_rows == 3);
    LOG(DEBUG) << "Generating mesh for " << kpt;
    km.k = kpt;
    km.kmax = kmax; // TODO: serves no purpose so far
    km.num = km.g.n_cols;
    km.r_K.set_size(km.num);
    // precalculate norm of k + G for all G
    for (size_t i = 0; i < km.num; ++i) {
        km.r_K(i) = arma::norm(bravais_mat*(km.k + km.g.col(i)));
    }
    LOG(DEBUG) << "r_K = \n" << print_mat(km.r_K);
    return km;
}

/** \brief Take the lower triangular matrices txx from Fleur and prepare them
 */
TMat::TMat(arma::cx_mat taa, arma::cx_mat tbb, arma::cx_mat tab, arma::cx_mat tba, int lnonsph)
{
    n = taa.n_rows;
     // each slice so far only contains the lower triangle of what are in total
    // 3 independent matrices:
    // AA and BB only need to be expanded
    aa = arma::symmatl(taa);
    bb = arma::symmatl(tbb);
    // AB and BA need to be combined in a single matrix
    arma::cx_vec diag_AB = tab.diag();
    ab =   arma::trimatl(tab)
                    + arma::trimatu(tba.t());
    ab.diag() = diag_AB;
    // Sane input data: Diagonals of inputs tab and tba must match
    arma::vec diff = arma::abs(diag_AB-tba.diag());
    bool diag_eq = arma::all(diff < 1e-15);
    assert(diag_eq);
    // ab and ba are really the same matrix, only hermitian transposed
    ba = ab.t();
    // calculate dimension from cutoff value
    L_nonspherical = (lnonsph+1)*(lnonsph+1);
    L_diag = n - L_nonspherical;
}

/** \brief Explicitly set T-matrices to zero for l > lmax_nonsph
 *
 * Rationale: The nonspherical (i.e. dense) part of the T-matrices is smaller
 * than the spherical (i.e. diagonal) part due to the different l-cutoffs.
 * Nevertheless, the T-matrices from Fleur contain the nonspherical parts up to the
 * larger spherical cutoff, so we explicitly zero them.
 */
void TMat::zero_nonsph(const int lmax_nonsph)
{
    int lm_start = (lmax_nonsph+1)*(lmax_nonsph+1);
    if (lm_start == n) {
        LOG(DEBUG) << "Nothing to zero, lsph = lnonsph.";
        return;
    }
    auto apply = [&](arma::cx_mat& T) {
        auto z = arma::zeros<arma::cx_mat>(n-lm_start, n);
        T.submat(lm_start, 0, n-1, n-1) = z;
        T.submat(0, lm_start, n-1, n-1) = z.t();
    };
    apply(aa);
    apply(bb);
    apply(ab);
    apply(ba);
}

/** \brief Overload for ostreams like cout. Only output T_aa, debug purposes
 */
std::ostream& operator<<(std::ostream &os, const TMat& tm)
{
    os << "T_aa =\n" << tm.aa << "\n";
       //<< "T_aa = " << tm.aa << "\n"
       //<< tm.bb << tm.ab << tm.ba;
    return os;
}

} // namespace flapw
