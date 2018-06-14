/** \file
 * Explicitly specialized template functions for different Binmat data types.
 */

#include "binmat.h"
namespace flapw {
template <>
bool Binmat<double>::check_type(const char c) { return c == 'R'; }
template <>
bool Binmat<std::complex<double>>::check_type(const char c) { return c == 'C'; }
template <>
bool Binmat<int>::check_type(const char c) { return c == 'I'; }

} // namespace flapw
