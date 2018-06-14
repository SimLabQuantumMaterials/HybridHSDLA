/** \file
 * Test the Fortran interop layer in libflerp.
 * TODO: comparison data for the meaningful tests
 */
#include "fleur_interop.h"
#include "binmat.h"
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <memory>


#include "test_common.h"

using namespace flapw;

/**
 * Test if the library works at all
 */
TEST(Basic)
{
    int a = 10;
    testfunc_(&a);
    CHECK(a == 30);
}


/**
 * Test if calling the wrapped sphbes function
 * works.
 */
TEST(Bessel)
{
    int lmax = 3;
    double x = 3.41;
    size_t length = lmax+1;
    //unique_ptr<double[]> jl(new double[30]);
    //double jl[length], djl[length];
    vector<double> jl(length, 0.), djl(length, 0.);
    wrap_bessel_(lmax, x, &jl[0], &djl[0]);
    cout << "jl = ";
    for (auto i = 0; i < length; ++i) 
        cout << jl[i] << " ";
    cout << endl;
    cout << "djl = ";
    for (auto i = 0; i < length; ++i) 
        cout << djl[i] << " ";
    cout << endl;
    CHECK(true);
}
/**
 * Test the wrapped ylm4 function from Fleur.
 */
TEST(SphericalYlm)
{
    int lmax = 3;
    double r[] = { 3.12, 4.3, 9.3 };
    size_t length = (lmax+1)*(lmax+1);
    vector<complex<double> > ylm(length);
    wrap_ylm_(lmax, r, &ylm[0]);
    cout << "ylm = ";
    for (auto i = 0; i < length; ++i)
        cout << ylm[i] << " ";
    cout << endl;
    CHECK(true);
}
/**
 * Test the wrapped radfun from Fleur. Relies on 
 * a specific dumped potential file for the hard-coded parameters.
 */
TEST(Radfun)
{
    string path = SI_POT_L0_PATH;
    int l = 0;
    double E = -0.145868090423439; 
    double dx = 2.2e-2;
    double r_MT = 2.16;
    RealBinmat vr;
    vr.read(path);
    int n_radpts = vr.get_size().first;
    auto sol = unique_ptr<double[]>(new double[2*n_radpts]);
    auto sol_dE = unique_ptr<double[]>(new double[2*n_radpts]);
    //double sol[n_radpts]
    //double sol_dE[19];
    double u, u_dr, dE_u, dE_u_dr, dE_u_norm;
    wrap_radfun_(l, E, dx, r_MT, &vr.get_data()[0], n_radpts, 
            sol.get(), sol_dE.get(),
            u, u_dr, dE_u, dE_u_dr, dE_u_norm);
    for (size_t i = 0; i < n_radpts; ++i) {
        cout << sol[i] << " ";
    }
    cout << endl;
    CHECK(true);
}

int main()
{
    return UnitTest::RunAllTests();
}
