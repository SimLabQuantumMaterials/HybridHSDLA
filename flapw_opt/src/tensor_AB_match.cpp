/** \file
 * Calculate the A- and B-tensors of matching coefficients
 */


#include "core.h"
#include "constants.h"
#include "special_funcs.h"
#include "types.h"
#include "AtomData.h"
#include <armadillo>

#include "logging.h"

namespace flapw {


/** \brief Generate the A and B tensors of matching coefficients.
 *
 * The formula: Same prefactor {1} for both A and B, indices {mu,L,G} suppressed.
 *
 * A = exp(iKp) 4pi/sqrt(\Omega ) 1/W i^l Y_{lm} (R_\mu K)    {1}
 *     [dE_u(rmt_a;l) |K| dj_l(rmt_a |K|) - dE_u_dr(rmt_a;l) j_l(rmt_a |K|)] {2}
 *
 * B = {1}
 *     [  -u(rmt_a;l) |K| dj_l(rmt_a |K|) +    u_dr(rmt_a;l) j_l(rmt_a |K|)] {3}
 *
 *  W = dE_u(rmt_a;l) u_dr(rmt_a;l) - u(rmt_a;l) dE_u_dr(rmt_a;l)
 *
 *
 * A and B are tensors with 3 dimensions (indices):
 *
 * * mu (atom index): 1...N_a (~100), number of atoms
 * * L  (sph. index): 1...(lmax+1)^2 (~9^2), enumeration of (l,m) tuples of spherical harmonics
 * * G  (rec. space): 1...(N_g*N_a) (~80*N_a), HUGE, enumerates all points in reciprocal space
 *
 * and a parameter:
 *
 * * k  (k-point): 3-vector, see below
 *
 * The parts explained:
 *
 * * K = k + G: `k` is supplied as a parameter, G ranges over all (precomputed) values s.th
 *     |k+G| < r_kmax, where r_kmax is a simulation parameter. 3-vector
 * * p: the position of the MT (i.e. position of the atom) in the global coordinate frame. 3-vector
 * * W: the Wronskian, arises from the solution of the coupling conditions. scalar
 * * Y_lm*: spherical harmonic for index L = (l,m), complex conjugate. scalar function of 3-vector
 * * R_mu: rotation matrix that transfers MT mu into the representative. 3x3 matrix
 * * u: solution of the radial equation for given `l` (from L) and atom type `alpha` for the l=0 component
 *      of the potential. similar symbols:
 *     *    dE_u: energy derivative
 *     *    u_dr: radial derivative
 *     * dE_u_dr: energy derivative of radial derivative
 * * rmt_a: radius of the muffin tin for atoms of type `alpha`. the index mu differs for each atom, but
 *          for atoms mu_1 and mu_2 of the same type (i.e. Silicon), alpha(mu_1) = alpha(mu_2). scalar
 * * |K|: 2-norm of vector K
 * * j_l: spherical bessel function for given `l`. scalar function of scalar.
 * * dj_l: derivative of j_l. TODO: fishy: derivative wrt. the scalar input? investigate.
 * * sqrt(\Omega): \Omega is the volume of the unit cell spanned by the real basis vectors
 *
 * \param   bmat        Basis vectors of reciprocal space supplied as its rows
*/
void construct_AB(arma::cx_cube &A, arma::cx_cube &B, const k_mesh& km,
                const AtomData& ad)

{
    // powers of imaginary i, i^l == i^(l mod 4)
    // FIXME: currently unused, as included in the tlm from Fleur
    cx_double i_power[4] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    int max_lmax = ad.max_lmax;
    const int size_lm = (max_lmax+1)*(max_lmax+1);
    A.zeros(size_lm, ad.num, km.num);
    B.zeros(size_lm, ad.num, km.num);
    /* Precompute solutions to the radial equation.
    * They are needed for each atom's MT radius and for all l values.
    * Only the value at the MT boundary is needed.
    * TODO:
        * lmax might (will?) depend on atom type
        ** by calculating for the overall lmax, additional computation occurs, but should be negligible
        * energy E usually depends on the l used (band), and atom type
        ** partial fix: l-dependence included
        * V is the same for atoms of the same type
    */
    auto u_match = ad.get_u_match();
    auto inv_wronk = ad.get_inv_wronk();
    auto vol_prefac = 4*pi/sqrt(ad.vol_cell);
    LOG(DEBUG) << "Precomputation finished.";
    for (size_t i = 0; i < km.num; ++i) {
        auto g = km.g.col(i);
        auto K = km.k + g;
        // rows of B are basis vectors: transpose and multipy with coefficients
        auto K_norm = km.r_K(i);
        // FIXME: rotation matrix for equivalent atoms, don't use equivalent
        // atoms, but that needs the potential for every atom (and not only representative)
        //auto K_base = arma::normalise(bmat.t() * K); // see gk in hssphn
        for (size_t j = 0; j < ad.num; ++j) {
            arma::vec3 K_base = arma::normalise(ad.get_rotmat(j).t() * K); // see gkrot in hssphn
            auto p = ad.pos.col(j); // position of current atom
            cx_double phase = std::exp(
                                    cx_double(0., 2*pi*arma::dot(K, p))
                                ); // phase factor exp(iKp)
            /* Precompute spherical Bessel functions and derivatives
             * for all l and given atom & K value.  */
            arma::mat bessel = calc_bessel(ad.get_mesh(j).rmt * K_norm, max_lmax);
            arma::cx_vec ylm_conj = arma::conj(calc_ylm(K_base, max_lmax));
            // assemble prefactor
            // no size_t on purpose, comparing signed/unsigned is a bad idea.
            for (int l = 0; l <= ad.get_lmax(j).first; ++l) {
                auto A_u = u_match(l, j).dE_u   *K_norm*bessel(l, 1)
                         - u_match(l, j).dE_u_dr       *bessel(l, 0);
                auto B_u = -u_match(l, j).u   *K_norm*bessel(l, 1)
                         +  u_match(l, j).u_dr       *bessel(l, 0);
                // FIXME: i_power commented out, because the Fleur guys already
                // include it in the T-matrices I need to dump. for the overlap,
                // it doesn't matter, since it cancels out anyway
                //auto common = vol_prefac*i_power[l % 4]*phase*inv_wronk(l,j);
                auto common = vol_prefac*phase*inv_wronk(l,j);
                auto sqrt_dE_u_norm = sqrt(u_match(l, j).dE_u_norm);
                for (int m = -l; m <= l; ++m) {
                    auto lm = l*(l+1)+m; // combined index for l and m
                    A(lm, j, i) = ylm_conj(lm) * common * A_u;
                    B(lm, j, i) = ylm_conj(lm) * common * B_u;
                }
            }
        }
    }
    LOG(DEBUG) << "Finished creating AB matrices of size "
               << A.n_rows << "," << A.n_cols << "," << A.n_slices;
}

void construct_A(arma::cx_cube &A, const k_mesh& km, const AtomData& ad) {
    cx_double i_power[4] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1} };
    int max_lmax = ad.max_lmax;
    const int size_lm = (max_lmax + 1) * (max_lmax + 1);
    A.zeros(size_lm, ad.num, km.num);
    auto u_match = ad.get_u_match();
    auto inv_wronk = ad.get_inv_wronk();
    auto vol_prefac = 4 * pi / sqrt(ad.vol_cell);
    LOG(DEBUG) << "Precomputation finished.";
    for (size_t i = 0; i < km.num; ++i) {
        auto g = km.g.col(i);
        auto K = km.k + g;
        auto K_norm = km.r_K(i);
        for (size_t j = 0; j < ad.num; ++j) {
            arma::vec3 K_base = arma::normalise(ad.get_rotmat(j).t() * K);
            auto p = ad.pos.col(j);
            cx_double phase = std::exp(cx_double(0, 2 * pi * arma::dot(K, p)));
            arma::mat bessel = calc_bessel(ad.get_mesh(j).rmt * K_norm, max_lmax);
            arma::cx_vec ylm_conj = arma::conj(calc_ylm(K_base, max_lmax));
            for (int l = 0; l <= ad.get_lmax(j).first; ++l) {
                auto A_u = u_match(l, j).dE_u * K_norm * bessel(l, 1)
                         - u_match(l, j).dE_u_dr * bessel(l, 0);
                auto common = vol_prefac * phase * inv_wronk(l,j);
                auto sqrt_dE_u_norm = sqrt(u_match(l, j).dE_u_norm);
                for (int m = -l; m <= l; ++m) {
                    auto lm = l * (l + 1) + m;
                    A(lm, j, i) = ylm_conj(lm) * common * A_u;
                }
            }
        }
    }
    LOG(DEBUG) << "Finished creating A matrix of size "
               << A.n_rows << "," << A.n_cols << "," << A.n_slices;
}

void construct_B(arma::cx_cube &B, const k_mesh& km, const AtomData& ad) {
    cx_double i_power[4] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1} };
    int max_lmax = ad.max_lmax;
    const int size_lm = (max_lmax + 1) * (max_lmax + 1);
    B.zeros(size_lm, ad.num, km.num);
    auto u_match = ad.get_u_match();
    auto inv_wronk = ad.get_inv_wronk();
    auto vol_prefac = 4 * pi / sqrt(ad.vol_cell);
    LOG(DEBUG) << "Precomputation finished.";
    for (size_t i = 0; i < km.num; ++i) {
        auto g = km.g.col(i);
        auto K = km.k + g;
        auto K_norm = km.r_K(i);
        for (size_t j = 0; j < ad.num; ++j) {
            arma::vec3 K_base = arma::normalise(ad.get_rotmat(j).t() * K);
            auto p = ad.pos.col(j);
            cx_double phase = std::exp(cx_double(0, 2 * pi * arma::dot(K, p)));
            arma::mat bessel = calc_bessel(ad.get_mesh(j).rmt * K_norm, max_lmax);
            arma::cx_vec ylm_conj = arma::conj(calc_ylm(K_base, max_lmax));
            for (int l = 0; l <= ad.get_lmax(j).first; ++l) {
                auto B_u = -u_match(l, j).u * K_norm * bessel(l, 1)
                         +  u_match(l, j).u_dr * bessel(l, 0);
                auto common = vol_prefac * phase * inv_wronk(l,j);
                auto sqrt_dE_u_norm = sqrt(u_match(l, j).dE_u_norm);
                for (int m = -l; m <= l; ++m) {
                    auto lm = l * ( l + 1) + m;
                    B(lm, j, i) = ylm_conj(lm) * common * B_u;
                }
            }
        }
    }
    LOG(DEBUG) << "Finished creating B matrix of size "
               << B.n_rows << "," << B.n_cols << "," << B.n_slices;
}

} // namespace flapw
