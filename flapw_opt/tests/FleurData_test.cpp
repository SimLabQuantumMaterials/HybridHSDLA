#include <iostream>
#include <type_traits>
#include "FleurData.h"
#include "types.h"

#include "test_common.h"

TEST(AtomDataPopulate)
{
    flapw::FleurData fd(KPT_DUMP_DIR);
    auto ad = fd.populate_atomdata();
    bool pos_eq = arma::all(arma::all(ad.pos == fd.get_atom_pos()));
    CHECK(pos_eq);
}
TEST(AtomData)
{
    flapw::AtomData ad;
    arma::mat p_tmp(3,5);
    p_tmp.randu();
    ad.pos = p_tmp;
    cout << "p_tmp = \n" << p_tmp << "\nad.pos = \n" << ad.pos;
    bool eq = arma::all(arma::all(ad.pos == p_tmp));
    CHECK(eq);
}

int main()
{
    return UnitTest::RunAllTests();
}
