#pragma once
/** \file
 * Provide a cleaner interface to manage data dumped from Fleur.
 */

#include "types.h"
#include "AtomData.h"
#include <map>
#include <string>
#include <armadillo>
#include <hdf5.h>

namespace flapw {

/** \brief Contains functions to retrieve data dumped from Fleur
 */
class FleurData {
public:
    ~FleurData();
    FleurData(const FleurData&) = delete;
    FleurData& operator=(const FleurData&) = delete;
    FleurData(std::string); ///< initialize the data directory

    std::pair<int,int> get_lmax(size_t i) const;
    arma::mat get_atom_pos() const;
    rad_mesh get_mesh(size_t type) const;
    arma::vec get_pot(size_t type) const;
    arma::vec get_energies(size_t type) const;
    arma::mat33 get_real_basis() const;
    arma::mat33 get_reciprocal_basis() const;
    arma::vec3 get_kpoint(size_t i) const;
    arma::imat get_g_mesh(size_t i) const;
    TMat get_fleur_tmats(size_t type) const;
    arma::mat33 get_rotmat(size_t atom) const;
    AtomData populate_atomdata(); ///< TODO: hackish way to fill the atomdata struct
    int get_num_kpts() const; 
protected:
    void read_scalars(const std::string p);
    void read_atom_pos(const std::string p);
    void read_pot_files(const std::string prefix);
    void read_kpt_dump(const std::string, const std::string);
    //void read_tmat(const std::string);
    void read_rotmat(const std::string);
    void hdf_atoms(const std::string);
    void hdf_basis(const std::string amat, const std::string bmat);
    void hdf_radmesh(const std::string);
    void hdf_lmax(const std::string loc_l, const std::string loc_lnonsph);
    void hdf_energies(const std::string);
    void hdf_potential(const std::string);
    void hdf_kpts(const std::string);
    void hdf_tmats(const std::string);
    int type_ind(size_t atom);
private:
    static const std::map<std::string, std::string> partmap; ///< Filename mappings
    // simulation data follows
    std::string data_folder; ///< root folder containing data files
    hid_t h5file; ///< file handle to the hdf data
    // "scalars"
    arma::ivec lmax; ///< maximum value of l, 1 entry per type
    arma::ivec lmax_nonsph; ///< like lmax, but for the nonspherical cutoff
    int num_atoms; ///< total number of individual atoms
    int num_types; ///< number of different atom types, according to fleur
    arma::ivec types; ///< num_atoms entries with the type of each atom
    arma::mat energies; ///< (lmax+1) x num_types energies. lmax is the maximum over all types' lmax.
    // vector, since it needs to be stored contiguously for HDF reading to work
    std::vector<rad_mesh> mesh_params; ///< mesh data for each atom's radial solution
    int max_mesh_pts; ///< maximum over all mesh_params' grid sizes
    int num_kpts; ///< total number of k-points
    int max_gpoints; ///< maximum number of g-points for any given k-vector
    arma::ivec num_gpts; ///< number of g-points for a given k-vector
    // "vectors"
    arma::mat positions; ///< each column is one atom's position
    arma::mat potential; ///< each col. is the potential for a value of l
    arma::mat33 real_basis; ///< column-wise basis vectors of real space
    arma::mat33 rec_basis; ///< row-wise basis vectors of reciprocal space
    arma::mat k_points; ///< columns are k-points and associated weight
    /// different # of g-pts for every k-point; g-vectors are cols. of the matrices
    arma::icube g_mesh;
    arma::field<TMat> t_matrices; ///< one t-matrix per slice
    arma::cube rotmats; ///< FIXME: rotation matrices for equivalent atoms
};

} // namespace flapw
