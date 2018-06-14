#ifndef FLEUR_INTEROP_H
#define FLEUR_INTEROP_H
/** \file
 * Definitions for functions wrapped from the
 * Fleur code base.
 * `const`-qualifiers mostly to distinguish what
 * is out- and input argument.
 */

#include <complex>

extern "C" {
void testfunc_(int*);
/** \brief Calculate the spherical bessel function & derivative
 *
 * Use the fleur function
 * \param       lmax    Value of l up and including to function values are returned
 * \param       x       Scalar argument to the Bessel function
 * \param[out]  jl      Size lmax+1 array with function values for l = 0..lmax
 * \param[out]  djl     Size lmax+1 array with *derivative* values
*/
void wrap_bessel_(
        const int& lmax, 
        const double& x, 
        double *const jl, 
        double *const djl
        ); 
/** \brief Calculate spherical harmonics
 *
 * Will return a flat array of all \f$Y_{lm}\f$ for \f$|m| \leq l,
 * l = 0\ldots \text{lmax}\f$, laid out in consecutive blocks for each l.
 * The blocks have size \f$2l+1\f$ and are thus
 *
 *     l = |  0  |    1   |    l = 2    | ....
 *     m = |  0  | -1 0 1 | -2 -1 0 1 2 | ....
 *
 * The value of \f$ Y_{lm}\f$ is found at index \f$ l(l+1)+m \f$.
 * The total array size is \f$(\text{lmax}+1)^2\f$, because (\f$ n:=\text{lmax}\f$)
 * \f[
 *  \sum_{l=0}^n 2l+1 = 2\cdot n(n+1)/2 + (n+1) = (n+1)(n+1) = (n+1)^2
 * \f]
 *
 * \param       lmax    Maximum value of l to generate function values for
 * \param       r       Vector of size 3, will be normalized and taken as input value
 * \param[out]  ylm     Complex, size (lmax+1)^2, see description
 */
void wrap_ylm_(
        const int& lmax,
        const double *const r, 
        std::complex<double> *const ylm
        );
/** \brief Solve the relativistic radial equation, produce matching coefficients
 *
 * The scalar relativistic approximation is used for the solution of the radial
 * Dirac equation on an exponential mesh.
 * Since the relativistic solution has two components, the respective arrays
 * have two columns (and are stored in column-major layout).
 *
 * Input values and basic solution:
 *
 * \param       l           radial parameter of the solution
 * \param       E           solve at that energy
 * \param       dx          grid spacing, i.e. \f$r_{i+1}=r_i \exp(\mathrm dx)\f$
 * \param       rmt         radius of outermost grid point, "muffin tin" radius
 * \param       vr          radial potential at each grid point
 * \param       n_radpts    number of gridpoints, must match size of vr
 * \param[out]  sol         large and small component of solution, n_radpts x 2
 * \param[out]  sol_dE      energy derivative of solution, same dimension
 *
 * Furthermore, special function values are returned separately:
 * The matching coefficients, all are real scalars and values of the "large"
 * component at the outermost grid point, the muffin tin (MT) radius "rmt".
 * They are needed to match the function inside the MT "sphere" C1-continuous
 * (i.e. function and derivative) to the outside.
 *
 * \param[out]  u           function value
 * \param[out]  u_dr        function's radial derivative
 * \param[out]  dE_u        energy derivative
 * \param[out]  dE_u_dr     radial derivative of energy derivative
 * \param[out]  dE_u_norm   Integral of energy derivative over whole mesh
 */
void wrap_radfun_(
        const int& l, const double& E, const double& dx, const double& rmt,
        const double *const vr, const int& n_radpts,
        double *const sol, double *const sol_dE,
        double& u, double& u_dr, double& dE_u, double& dE_u_dr,
        double& dE_u_norm
        );
}


#endif 
