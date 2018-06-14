#include "binmat.h"
#include "easylogging++.h"
#include <vector>
#include <string>

#include "test_common.h"

using namespace flapw;

vector<string> paths;

TEST(Readmatrix)
{
    for (auto p : paths) {
        cout << "Reading file " << p << ":\n";
        char c = RealBinmat::file_test_type(p);
        IntBinmat ibm;
        RealBinmat rbm;
        ComplexBinmat cbm;
        switch(c) {
            case 'I':
                ibm.read(p);
                cout << ibm.get_data();
                break;
            case 'R':
                rbm.read(p);
                cout << rbm.get_data();
                break;
            case 'C':
                cbm.read(p);
                cout << cbm.get_data();
                break;
        }
    }
    cout << "END" << endl;
    CHECK(true);
}

int main(int argc, char** argv) {
    if (argc < 2) { // standard testing mode
        paths.push_back(SI_POT_L0_PATH);
        paths.push_back(BINMAT_GPT_PATH);
    } else { // user-supplied matrix
        paths.push_back(argv[1]);
    }
    return UnitTest::RunAllTests();
}
