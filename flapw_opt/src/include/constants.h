#pragma once
/** \file 
 * Includes definitions of physical constants.
 */

namespace flapw {
    // from CODATA 2010 web page, in SI units

    const double c_light = 299792458; /**< speed of light, m/s */
    const double el_mass = 9.10938291e-31; /**< electron mass, kg */
    const double el_charge = 1.602176565e-19; /**< electron (elementary) charge */
    const double fs_alpha = 137.035999074; /**< 1/(fine-structure constant alpha), dimensionless */
    const double fs_alpha_fleur = 137.0359895; /**< as above, but "old" value from fleur */

    // math constants
    const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

    /** Different unit systems produce different values for some constants */
    enum unit_system { HARTREE, RYDBERG };

    const double c_atomic[] = { [HARTREE] = fs_alpha_fleur, [RYDBERG] = 2*fs_alpha_fleur };

}
