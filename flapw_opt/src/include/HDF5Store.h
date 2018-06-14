#pragma once
/** \file
 * Provide convenience methods for storing simulation results into a HDF5 file
 */

#include <string>
#include <hdf5.h>
#include <armadillo>


class HDF5Store {
public:
    HDF5Store(const std::string filename);
    ~HDF5Store();
    void add_hsmat(const arma::cx_mat& H, const arma::cx_mat& S, const int k_index);
    void add_arma(const arma::cx_cube& M, const std::string& name);
    void add_arma(const arma::cx_mat& M, const std::string& name);
    void flush();
protected:
    void write_cx_data(int rank, const hsize_t dims[], const arma::cx_double* pM, const std::string& name);
private:
    hid_t h5file;
    hid_t g_hmat;
    hid_t g_smat;
    hid_t h5cx;
};
