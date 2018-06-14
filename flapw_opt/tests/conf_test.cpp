#include "constants.h"
#include <iostream>

#include "test_common.h"

using namespace flapw;

TEST(Constants)
{
    double c = c_light;
    cout << c << endl;
    double c_ry = c_atomic[RYDBERG];
    double c_ht = c_atomic[HARTREE];
    CHECK(c_ry == 2*c_ht);
}

int main() {
    return UnitTest::RunAllTests();
}
