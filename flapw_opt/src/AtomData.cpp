#include "AtomData.h"
#include "special_funcs.h"
#include "logging.h"
#include "util.h"
#include <cassert>

namespace flapw {

/** \brief TBD
 */
AtomData::AtomData()
: tmats_augmented(false), exists_radial_solution(false)
{

}

/** \brief Prepare the B-tensor to use in Overlap matrix generation
 *
 * This makes a full copy of B_no_norm (i.e. the "standard" B) and applies the
 * factor of the matching coefficient as a square root, so the whole operation can
 * still be expressed as a HERK
 */

void AtomData::apply_u_norm(arma::cx_cube& B) const
{
    for (int at = 0; at < num; ++at) {
        int lmax = get_lmax(at).first;
        for (int l = 0; l <= lmax; ++l) {
            auto sqrt_dE_u_norm = sqrt(u_match(l, at).dE_u_norm);
            for (int m = -l; m <= l; ++m) {
                int lm = l*(l+1) + m;
                B.tube(lm, at) *= sqrt_dE_u_norm;
            }
        }
    }
}


/** \brief Return mesh data for a given atom */
rad_mesh AtomData::get_mesh(size_t atom) const
{
    return mesh[atom];
}

/** \brief Return the potential of an atom */
arma::vec AtomData::get_pot(size_t atom) const
{
    return V_rad[atom];
}


/** \brief Return l-cutoffs for an atom
 *
 * \returns a pair, first entry is the spherical, second entry
 *          the nonspherical cutoff
 */
AtomData::l_pair AtomData::get_lmax(size_t atom) const
{
    return lmax[atom];
}

/** \brief Calculate the inverse Wronskian from given matching data */
inline double calc_inv_wronskian(const rad_match& c)
{
    return 1./(c.dE_u*c.u_dr - c.u*c.dE_u_dr);
}

/** \brief Generate matching coefficients and Wronskian
 *
 * Both quantities are needed for the creation of the AB-tensors
 */
void AtomData::generate_matchdata()
{
    if (exists_radial_solution) return; // does not need to be recalculated
    u_match.set_size(max_lmax+1, num); // FIXME: cannot easily set to 0, memset?
    inv_wronk.zeros(max_lmax+1, num); // TODO: this should (hopefully) prevent the "extra"
                                      // entries for atoms with smaller lmax from doing any damage

    for (size_t i = 0; i < num; ++i) {
        const auto E = get_E(i);
        for (int l = 0; l <= get_lmax(i).first; ++l) {
            radial_solve(l, E(l), get_mesh(i), get_pot(i), u_match(l, i));
            inv_wronk(l, i) = calc_inv_wronskian(u_match(l,i));
        }
    }
    LOG(DEBUG) << "Inverse Wronskian:\n" << inv_wronk;
    exists_radial_solution = true;
}

/** \brief Add the diagonal "spherical" terms
 *
 * They are considered separately in Fleur and not contained in the T-matrices
 * that are dumped.
 */
void AtomData::augment_tmats()
{
    // FIXME: really bad design, should move all the logic into a constructor, but no time
    // at least there is some error checking!
    if (t_matrices.size() == 0) {
        LOG(ERROR) << "Cannot augment T-matrices, none are loaded";
        return;
    }
    if (tmats_augmented) {
        LOG(DEBUG) << "T-Matrices are already augmented";
        return;
    }
    generate_matchdata();
    assert(t_matrices.size() == num);
    /* FIXME or HACK, you decide...
     * Fleur has the distinction between spherical and nonspherical terms...
     * It does: 1) loop only over diagonal up to lmax, apply diagonal terms (cheap); those are added here to T
     * 2) loop over full matrix, but only to lmax_nonsph < lmax, apply terms contained in T
     * In principle, using a single cutoff should be more accurate, and with the method implemented here
     * (just a GEMM in terms of matrices, no separation as in Fleur) it is easier to implement.
     * Though for consistency with Fleur, we replicate its behavior (even if it might not make sense numerically).
     * Better solutions (?):
     * 1) explicitly treat spherical terms separately
     * 2) use single, larger cutoff - sacrifice consistency with Fleur since we'll be "better"
     */
    for (size_t atom = 0; atom < num; ++atom) {
        auto& T = t_matrices[atom];
        const auto& El = energies[atom];
        const auto& rm = u_match.col(atom);
        // zero T for nonspherical_l+1 ... spherical l cutoff so final result only has diagonal
        T.zero_nonsph(lmax[atom].second);
        // 0...spherical l cutoff
        for (int l = 0; l <= lmax[atom].first; ++l) {
            for (int m = -l; m <= l; ++m) {
                int lm = l*(l+1) + m;
                T.aa(lm, lm) += El(l);
                T.bb(lm, lm) += El(l)*rm(l).dE_u_norm;
                // do as Fleur does, don't do as Fleur says...
                // symmetrize BA and AB explicitly
                T.ab(lm, lm) += .5;
                T.ba(lm, lm) += .5;
            }
        }
        assert(is_hermitian(T.aa));
        assert(is_hermitian(T.bb));
        bool AB_is_BA_trans = arma::all(arma::all(T.ab == T.ba.t()));
        assert(AB_is_BA_trans);
    }
    tmats_augmented = true;
    LOG(DEBUG) << "Done augmenting T-matrices with diagonal terms.";
}

/** \brief Return a reference to the internal storage of all T-matrices
 *
 * Returns a const ref. to avoid copying all the matrices.
 * CAREFUL: The reference will go dangling if AtomData goes out of scope.
 */
const AtomData::TMatVec& AtomData::get_tmats_handle() const
{
    return t_matrices;
}


match_field AtomData::get_u_match() const
{
    if (!exists_radial_solution) throw "Need to call generate_matchdata() first.";
    return u_match;
}
arma::mat AtomData::get_inv_wronk() const
{
    if (!exists_radial_solution) throw "Need to call generate_matchdata() first.";
    return inv_wronk;
}


/** \brief Return Energy parameters for a given atom */
arma::vec AtomData::get_E(size_t atom) const
{
    return energies[atom];
}

arma::mat33 AtomData::get_rotmat(size_t atom) const
{
    return rotmats[atom];
}



}
