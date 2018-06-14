#pragma once
/** \file
 * Again, for lack of a better name: Core function as in "this is what
 * we really want to calculate, i.e. the AB-tensors for now
 */

#include "types.h"
#include "AtomData.h"
#include <armadillo>

using arma::cx_double;

namespace flapw {

void construct_AB(arma::cx_cube &A, arma::cx_cube &B, const k_mesh& km,
                const AtomData& ad);

void construct_A(arma::cx_cube &A, const k_mesh& km, const AtomData& ad);

void construct_B(arma::cx_cube &B, const k_mesh& km, const AtomData& ad);

} // namespace flapw
