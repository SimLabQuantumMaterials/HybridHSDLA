#pragma once

#include "types.h"
#include <armadillo>

namespace flapw {

arma::mat radial_solve(const int l, const double E, const rad_mesh& rm,
            const arma::vec& V, rad_match& coeff);
arma::mat calc_bessel(const double x, const int lmax);
arma::cx_vec calc_ylm(const arma::vec3& r, const int lmax);

} // namespace flapw
