#pragma once
/** \file
 * AtomData is the interface class storing atom-specific information such as
 * meshes, potential, ...
 * It is agnostic of equivalent atoms and stores the data for each atom separately.
 */

#include "types.h"

namespace flapw {

/** Store information about atoms in the simulation cell
 *
 *  Positions, types, (symmetries?), ...
 */
class AtomData {
public:
    using l_pair = std::pair<int, int>;
    using TMatVec = std::vector<TMat>;
    AtomData();

    size_t num; ///< number of atoms
    int max_mesh_pts; ///< maximum value of mesh[i].num
    int max_lmax; ///< maximum spherical l-cutoff
    double vol_cell; ///< volume of the cell
    arma::mat pos; ///< position of atom i in column i
    std::vector<rad_mesh> mesh; ///< mesh information, 1 entry per atom
    std::vector<arma::vec> V_rad; ///< potential data, 1 entry per atom
    std::vector<l_pair> lmax; ///< spherical and nonspherical l-cutoffs
    std::vector<arma::vec> energies; ///< per-atom energy parameters
    std::vector<arma::mat33> rotmats; ///< per-atom rotation matrices
    std::vector<TMat> t_matrices; ///< complete T-matrices according to the formulae in pkurz
    // TODO
    // radial solution storage
    match_field u_match;
    arma::mat inv_wronk;
    void generate_matchdata();
    void augment_tmats();
    void apply_u_norm(arma::cx_cube& B) const;
    match_field get_u_match() const;
    arma::mat get_inv_wronk() const;
    // ?type
    rad_mesh get_mesh(size_t atom) const;
    arma::vec get_pot(size_t atom) const;
    auto get_lmax(size_t atom) const -> l_pair;
    arma::vec get_E(size_t atom) const;
    arma::mat33 get_rotmat(size_t atom) const;
    const TMatVec& get_tmats_handle() const;
private:
    bool tmats_augmented;
    bool exists_radial_solution;
};

} // namespace flapw
