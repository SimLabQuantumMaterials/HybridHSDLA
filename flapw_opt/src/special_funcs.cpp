/** \file
 * For lack of a better name: "Special" mathematical functions that provide
 * data for the AB-tensor construction (spherical harmonics, ...)
 */

#include "special_funcs.h"
#include "types.h"
#include "fleur_interop.h"
#include <armadillo>
#include <cassert>

#include "logging.h"

namespace flapw {

/** Solve the radial equation
 *
 * TODO: is E an atom-specific parameter that should go into rad_mesh?
 * l definitely is not, it is specific to the individual solution.
 */
arma::mat radial_solve(const int l, const double E, const rad_mesh& rm,
            const arma::vec& V, rad_match& coeff)
{
    assert(V.n_elem == rm.num); // potential and solution mesh must match
    // i-th row is evaluated at i-th grid point
    // each column is a separate component of the (relativistic) solution:
    // 0   |  1   |  2  |  3
    // big | small| big | small
    // solution   | energy derivative
    arma::mat solution(rm.num, 4);
    wrap_radfun_(l, E, rm.dx, rm.rmt, V.memptr(), rm.num,
            solution.colptr(0), solution.colptr(2),
            coeff.u, coeff.u_dr, coeff.dE_u, coeff.dE_u_dr,
            coeff.dE_u_norm);
    return solution;
}



/** Calculate the sperical bessel function and its derivative.
 *
 * For a given value of lmax, lmax+1 values will be returned.
 * The first column contains j_l, the second one dj_l (the derivative).
 */
arma::mat calc_bessel(const double x, const int lmax)
{
    arma::mat res(lmax+1, 2);
    wrap_bessel_(lmax, x, res.colptr(0), res.colptr(1));
    return res;
}
/** \brief Calculate all spherical harmonics up to \f$ l = \text{lmax}\f$
 *
 * For a description of the layout of entries, \see wrap_ylm_()
 */
arma::cx_vec calc_ylm(const arma::vec3& r, const int lmax)
{
    arma::cx_vec res((lmax+1)*(lmax+1));
    wrap_ylm_(lmax, r.memptr(), res.memptr());
    return res;
}

} // namespace flapw
