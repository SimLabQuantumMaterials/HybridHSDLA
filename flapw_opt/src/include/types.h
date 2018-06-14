#pragma once
/** \file
 * Simulation data types and helper prototypes.
 */

#include <armadillo>

namespace flapw {

/** Store mesh of g-values around a specified k-point
 *
 * All G satisfy \f$|k+G| < k_\text{max}\f$, i.e. lie within a "sphere"
 */
struct k_mesh {
    double kmax; ///< "radius" of the sphere
    arma::vec3 k; ///< mesh generated around this k-point
    arma::ivec3 max_ind; ///< maximum index in each dimension
    arma::imat g; ///< hold all valid g-points, size 3 x num
    arma::vec r_K; ///< radius_K, i.e. norm of k + G forall G
    size_t num; ///< number of g-points, i.e. g.n_cols
};

/** Information to solve the radial equation on an exponential mesh.
 */
struct rad_mesh {
    int num; ///< number of mesh points
    double dx; ///< mesh increment: r_i+1 = r_i * exp(dx)
    double rmt; ///< muffin-tin radius (outermost mesh point)
};

/** Used to return matching coefficients from the radial solver
 *
 * Represent values of function/derivatives at the MT radius that was passed in.
 */
struct rad_match {
    double u;  ///< function itself
    double u_dr; ///< radial derivative
    double dE_u; ///< energy derivative function
    double dE_u_dr; ///< radial derivative of energy derivative function
    /** not a matching coefficient: norm (i.e. integral) of the energy derivative
     * function over the whole MT radial grid. */
    double dE_u_norm;
};

using match_field = arma::field<rad_match>; ///< matching data for (lmax+1, atoms)

/** \brief Storage for a combined set of 4 T-matrices */
struct TMat {
    TMat() : n(0) { };
    TMat(arma::cx_mat taa, arma::cx_mat tbb, arma::cx_mat tab, arma::cx_mat tba, int lnonsph);
    size_t n; ///< size of the matrices
    int L_nonspherical; ///< dimension of the dense, nonspherical submatrix
    int L_diag; ///< length of the spherical-only diagonal part
    arma::cx_mat aa, bb, ab, ba;
    friend std::ostream& operator<<(std::ostream &os, const TMat& tm);
    void zero_nonsph(const int lmax_nonsph);
};

enum TMatType {T_AA, T_BB, T_AB, T_BA}; ///< used as array indices to select t-matrix

/* Non-class function prototypes */
k_mesh create_kmesh(const arma::vec3& kpt, const double kmax,
                    const arma::mat33& bravais_mat, arma::imat g_dump);

} // namespace flapw
