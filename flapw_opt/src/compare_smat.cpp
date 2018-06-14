#include "binmat.h"
#include "easylogging++.h"
#include <vector>
#include <string>
#include <armadillo>
#include <cassert>
#include <cstdlib> // for EXIT_FAILURE

_INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace flapw;
using namespace arma;

template <typename T>
Mat<T> op_minus(const Mat<T>& a, const Mat<T>& b) { return a - b;}

int main(int argc, char** argv) {
    vector<string> args;
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " arg1 arg2 [output]" << endl
             << "calculates 'arg1 - arg2' and optionally saves to output." << endl;
        return EXIT_FAILURE;
    }
    // TODO: only for complex
    ComplexBinmat bm1(argv[1]), bm2(argv[2]);
    auto res_pack = op_minus(bm1.get_data(), bm2.get_data());
    // reshape packed storage
    int n = (1 + sqrt(1 + 8*res_pack.n_elem))/2;
    cout << "n: " << n << ", n_elem: " << res_pack.n_elem << endl;
    assert( n*(n-1)/2 == res_pack.n_elem ); // matrix needs to be square
    cx_mat res(n, n);
    size_t row = 0, col = 0;
    for (auto el : res_pack) {
        res(row, col) = el;
        if (row == col) {
            ++row;
            col = 0;
        } else {
            ++col;
        }
    }

    if (argc > 3) { // output to file, ".mat" is ASCII and ".bin" is binary
        res.save(string(argv[3]) + ".mat", arma_ascii);
        res.save(string(argv[3]) + ".bin", arma_binary);
    } else {
        cout << res;
    }    
    return 0;
}
