#pragma once
/** \file
 * Small helper functions, for...
 * printing
 */

#include <iostream>
#include <armadillo>

template<typename T>
std::string print_mat(T mat)
{
    std::ostringstream os;
    for (auto p : mat) {
        os << p << " ";
    }
    return os.str();
}

/** \brief Check if the given matrix is hermitian to machine precision
 */
template<typename T>
bool is_hermitian(const T& mat)
{
    return all(all(
        arma::abs(mat - arma::trans(mat)) < 1e-15));
}
