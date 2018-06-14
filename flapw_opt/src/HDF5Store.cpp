#include "HDF5Store.h"
#include "logging.h"

/** \brief Initialize store with given filename
 */
HDF5Store::HDF5Store(const std::string filename)
{
    h5file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    LOG(DEBUG) << "Opened file " << filename;
    g_hmat = H5Gcreate2(h5file, "/matrix_hamilton", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    g_smat = H5Gcreate2(h5file, "/matrix_overlap", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // need to create our own complex type
    hsize_t dims[] = {2};
    h5cx = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
}

/** \brief Add a new Hamilton and Overlap matrix for the specified k-point index
 */
void HDF5Store::add_hsmat(const arma::cx_mat& H,
                          const arma::cx_mat& S,
                          const int k_index)
{
    bool square = H.is_square() && S.is_square();
    auto n = H.n_rows;
    if ((n != S.n_rows) || !square) {
        LOG(ERROR) << "H and S have different sizes / are not square: n = " << n;
    }
    hsize_t dims[] = {n, n};
    hid_t dsp = H5Screate_simple(2, dims, NULL);
    // add 1 so that indices are the same as in the fleur dump file
    auto kind_str = std::to_string(k_index + 1);
    {
        // H-matrix
        hid_t dset = H5Dcreate2(g_hmat, kind_str.c_str(), h5cx, dsp,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, h5cx, H5S_ALL, H5S_ALL, H5P_DEFAULT, H.memptr());
        H5Dclose(dset);
    }
    {
        // S-matrix
        hid_t dset = H5Dcreate2(g_smat, kind_str.c_str(), h5cx, dsp,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, h5cx, H5S_ALL, H5S_ALL, H5P_DEFAULT, S.memptr());
        H5Dclose(dset);
    }
    H5Sclose(dsp);
}

void HDF5Store::flush()
{
    H5Fflush(h5file, H5F_SCOPE_GLOBAL);
}


/** \brief Save a complex matrix at the root of the file
 */
void HDF5Store::add_arma(const arma::cx_cube& M, const std::string& name)
{
    int rank = 3;
    hsize_t dims[] = {M.n_rows, M.n_cols, M.n_slices};
    write_cx_data(rank, dims, M.memptr(), name);
}
void HDF5Store::add_arma(const arma::cx_mat& M, const std::string& name)
{
    int rank = 2;
    hsize_t dims[] = {M.n_rows, M.n_cols};
    write_cx_data(rank, dims, M.memptr(), name);
}

void HDF5Store::write_cx_data(int rank, const hsize_t dims[],
                              const arma::cx_double* pM, const std::string& name)
{
    hid_t dsp = H5Screate_simple(rank, dims, NULL);
    hid_t dset = H5Dcreate2(h5file, name.c_str(), h5cx, dsp,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, h5cx, H5S_ALL, H5S_ALL, H5P_DEFAULT, pM);
    H5Dclose(dset);
    H5Sclose(dsp);
}

HDF5Store::~HDF5Store()
{
    H5Tclose(h5cx);
    H5Gclose(g_hmat);
    H5Gclose(g_smat);
    H5Fclose(h5file);
}
