#include "FleurData.h"
#include "binmat.h"
#include "AtomData.h"
#include <fstream>
#include <cassert>
#include <exception>
// for zero-padding
#include <sstream>
#include <iomanip>
// new data format
#include <hdf5.h>
#include <hdf5_hl.h>

#include "easylogging++.h"

namespace flapw {

/**
 * Hard-coded prefixes/names of data files to be loaded from a given data folder
 * during object creation.
 */
const std::map<std::string, std::string> FleurData::partmap = {
    {"scalars",         "scalars.dat"},
    {"pot_prefix",      "radpot_l-"},
    {"k-points",        "kpts"},
    {"gpoint_prefix",   "gpts_"},
    {"positions",       "positions.bin"},
    {"amat",            "real_basis.bin"},
    {"bmat",            "rec_basis.bin"},
    {"tmat_prefix",     "tmat_"},
    {"rotmat_prefix",   "rotmat_"},
    {"hdf_name",        "fleur_dump.h5"}
};
/** \brief Return stored data in AtomData format
 *
 * Use the loaded data to construct an AtomData object that can be passed to
 * other functions.
 * \returns Fully setup AtomData object with this instance's data
 */
AtomData FleurData::populate_atomdata()
{
    AtomData ad;
    ad.num = num_atoms;
    ad.max_mesh_pts = max_mesh_pts;
    ad.pos = positions;
    ad.vol_cell = arma::det(real_basis);
    ad.max_lmax = lmax.max();
    assert( (potential.n_rows == max_mesh_pts) && (potential.n_cols == num_types) );
    // Replicate "raw" data and create vectors with one entry per atom
    // This discards information about equivalent atoms and instead stores it
    // for each atom individually
    for (size_t atom = 0; atom < num_atoms; ++atom) {
        int at_type = type_ind(atom);
        LOG(DEBUG) << "Atom " << atom << " has type " << at_type;
        //LOG(INFO) << "Atom " << atom << " has type " << at_type;
        ad.V_rad.push_back(get_pot(at_type));
        ad.mesh.push_back(get_mesh(at_type));
        ad.lmax.push_back(get_lmax(at_type));
        ad.energies.push_back(get_energies(at_type));
        ad.rotmats.push_back(get_rotmat(atom));
        ad.t_matrices.push_back(get_fleur_tmats(at_type));
    }
    assert( (ad.V_rad.size() == ad.num) );
    // augment the T-matrices: right now, they are missing the diagonal terms
    ad.augment_tmats();
    return ad;
}
/** \brief Translate atom index into a type index
 *
 * Types start from 1 in Fleur, but we need a 0-based array index
 */
int FleurData::type_ind(size_t atom)
{
    return types(atom) - 1;
}


/** \brief max value of l for radial solves/summation
 */
std::pair<int,int> FleurData::get_lmax(size_t i) const
{
    return std::make_pair(lmax(i), lmax_nonsph(i));
}
/** \brief Return the matrix of atom positions
 *
 * 3xN, each column contains the xyz coordinates of an atom
 */
arma::mat FleurData::get_atom_pos() const
{
    return positions;
}

/** \brief Meturn mesh data for given atom type
 */
rad_mesh FleurData::get_mesh(size_t type) const
{
    return mesh_params[type];
}

/** \brief Potential for given atom type
 */
arma::vec FleurData::get_pot(size_t type) const
{
    int num_mesh_pts = get_mesh(type).num;
    return potential(arma::span(0, num_mesh_pts-1), type);
}

/** \brief Return energies for all l for a given type
 */
arma::vec FleurData::get_energies(size_t type) const
{
    int lmax_spherical = get_lmax(type).first;
    return energies(arma::span(0, lmax_spherical), type);
}

/** \brief Return the basis matrix of real space
 *
 * Columns are basis vectors
 */
arma::mat33 FleurData::get_real_basis() const
{
    return real_basis;
}

/** \brief Return the basis of reciprocal space
 *
 * rows are basis vectors of rec. space
 */
arma::mat33 FleurData::get_reciprocal_basis() const
{
    return rec_basis;
}

/** \brief Return i-th k-point
 */
arma::vec3 FleurData::get_kpoint(size_t i) const
{
    return k_points.col(i);
}

/** \brief Return g-mesh for given k-point index
 */
arma::imat FleurData::get_g_mesh(size_t i) const
{
    std::cout << "get_g_mesh: " << i << "  " << num_gpts(i) << std::endl;
    using arma::span;
    return g_mesh(span(), span(0, num_gpts(i)-1), span(i));
}

/** \brief Return fleur's T-matrices for a given atom type, see TMat
 */
TMat FleurData::get_fleur_tmats(size_t type) const
{
    return t_matrices(type);
}

/** \brief Return rotation matrices for all atoms
 */
arma::mat33 FleurData::get_rotmat(size_t atom) const
{
    return rotmats.slice(atom);
}

/** \brief Return number of stored k-pts
 */
int FleurData::get_num_kpts() const
{
    return num_kpts;
}


/*******************
 * SETUP FUNCTIONS *
 *******************/
FleurData::FleurData(std::string folder) :
data_folder(folder), num_types(-1),
max_mesh_pts(0), num_kpts(0), max_gpoints(0)
{
    using std::string;
    /// Join path 'a' and 'b' by the path separator
    auto j = [](const string& a, const string& b) {return a + "/" + b; };
    string h5path = j(folder, partmap.at("hdf_name"));
    h5file = H5Fopen(h5path.c_str(),
                           H5F_ACC_RDONLY, H5P_DEFAULT);
    if (sizeof(arma::sword) == sizeof(int)) {
	h5_int_mem_type = H5T_NATIVE_INT;
    } else if (sizeof(arma::sword) == sizeof(long)) {
	h5_int_mem_type = H5T_NATIVE_LONG;
    } else if (sizeof(arma::sword) == sizeof(long long)) {
	h5_int_mem_type = H5T_NATIVE_LLONG;
    } else {
	LOG(FATAL) << "Could not detect int datatype.";
    }
    if (h5file < 0) {
        LOG(FATAL) << "Could not open file: " << h5path;
    }
    LOG(INFO) << "Opened HDF5 data file at " << h5path;
    hdf_basis("/real_basis", "/rec_basis");
    hdf_atoms("/atoms");
    hdf_radmesh("/mesh");
    hdf_lmax("/lmax", "/lmax_nonspherical");
    hdf_energies("/energy_params");
    hdf_potential("/potential");
    hdf_kpts("/k_meshes");
    hdf_tmats("/tmat");
    /*read_scalars(j(folder, partmap.at("scalars")));
    read_atom_pos(j(folder, partmap.at("positions")));
    read_pot_files(j(folder, partmap.at("pot_prefix")));
    read_basis(j(folder, partmap.at("amat")),
                      j(folder, partmap.at("bmat")));
    read_kpt_dump(j(folder, partmap.at("k-points")),
                  j(folder, partmap.at("gpoint_prefix")));
    read_tmat(j(folder, partmap.at("tmat_prefix")));
    read_rotmat(j(folder, partmap.at("rotmat_prefix")));
    */
}
FleurData::~FleurData()
{
    H5Fclose(h5file);
}
/** \brief Read t-matrices for all atom types (4 each)
 */
void FleurData::hdf_tmats(const std::string group_loc)
{
    hid_t g_tmat = H5Gopen2(h5file, group_loc.c_str(), H5P_DEFAULT);
    t_matrices.set_size(num_types);
    // type dataset names start at 1
    for (int type = 1; type <= num_types; ++type) {
        auto type_str = std::to_string(type);
        int n = 0;
        hid_t dset = H5Dopen(g_tmat, type_str.c_str(), H5P_DEFAULT);
        hid_t dtype = H5Dget_type(dset);
        // read size of current tmat-set and initialize temp storage
        H5LTget_attribute_int(dset, ".", "n", &n);
        //cx_mat taa(n,n), tbb(n,n), tab(n,n), tba(n,n);
        arma::cx_cube t_all(n, n, 4);
        H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, t_all.memptr());
        // TODO: check ordering, here as well as in fleur. maybe add enum to h5file?
        t_matrices(type-1) = TMat(t_all.slice(T_AA), t_all.slice(T_BB),
                                  t_all.slice(T_AB), t_all.slice(T_BA),
                                  get_lmax(type-1).second
                                 );
        H5Tclose(dtype);
        H5Dclose(dset);
    }
    H5Gclose(g_tmat);
    LOG(DEBUG) << "T-matrices for atom type 1:\n"
               << "[see file tmat_hdf_dump_aa0.bin]";
    t_matrices(0).aa.save("tmat_hdf_dump_aa0.bin", arma::arma_binary);
}

/** \brief Read k-point data and associated meshes
 */
void FleurData::hdf_kpts(const std::string loc)
{
    // get size first to allocate space for the data
    H5LTget_attribute_int(h5file, loc.c_str(), "max_gpts", &max_gpoints);
    H5LTget_attribute_int(h5file, loc.c_str(), "num_kpts", &num_kpts);
    num_gpts.set_size(num_kpts); // one entry for each k-pt
    k_points.set_size(3, num_kpts);
    g_mesh.set_size(3, max_gpoints, num_kpts);
    hid_t dset = H5Dopen(h5file, loc.c_str(), H5P_DEFAULT);
    hid_t dtype = H5Dget_type(dset);
    // read fields individually
    {
        // ---- actual number of g-pts ----
        const char* fieldname = "num_gpts";
        hid_t t_memb = H5Tget_member_type(dtype, H5Tget_member_index(dtype, fieldname));
        //hid_t t_field = H5Tcreate(H5T_COMPOUND, H5Tget_size(t_memb));
	LOG(INFO) << "size: "<< H5Tget_size(t_memb);
        hid_t t_field = H5Tcreate(H5T_COMPOUND, sizeof(arma::sword)); // DIEGO
        H5Tinsert(t_field, fieldname, 0,  h5_int_mem_type);
        //H5Tinsert(t_field, fieldname, 0, t_memb);
        //H5Dread(dset, h5_int_mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, num_gpts.memptr());
        H5Dread(dset, t_field, H5S_ALL, H5S_ALL, H5P_DEFAULT, num_gpts.memptr());
        H5Tclose(t_memb);
        H5Tclose(t_field);
        LOG(DEBUG) << "# gpts:\n" << num_gpts;
        //LOG(INFO) << "# gpts:\n" << num_gpts; // DIEGO
    }
    {
        // ---- the k-point ----
        const char* fieldname = "k_point";
        hid_t t_memb = H5Tget_member_type(dtype, H5Tget_member_index(dtype, fieldname));
        hid_t t_field = H5Tcreate(H5T_COMPOUND, H5Tget_size(t_memb));
        H5Tinsert(t_field, fieldname, 0, t_memb);
        H5Dread(dset, t_field, H5S_ALL, H5S_ALL, H5P_DEFAULT, k_points.memptr());
        H5Tclose(t_memb);
        H5Tclose(t_field);
        LOG(DEBUG) << "k-points:\n" << k_points;
    }
    {
        // ---- the mesh of g-points associated with each k-point ----
        const char* fieldname = "g_mesh";
        hid_t t_memb = H5Tget_member_type(dtype, H5Tget_member_index(dtype, fieldname));
        hid_t t_field = H5Tcreate(H5T_COMPOUND, H5Tget_size(t_memb));
        H5Tinsert(t_field, fieldname, 0, t_memb);
        H5Dread(dset, t_field, H5S_ALL, H5S_ALL, H5P_DEFAULT, g_mesh.memptr());
        H5Tclose(t_memb);
        H5Tclose(t_field);
#ifndef NDEBUG
        // do minimal pretty-printing on the g-mesh of the first k-point,
        // lots of data to print
        std::ostringstream os;
        auto slice0 = g_mesh.slice(0);
        for (size_t g = 0; g < slice0.n_cols; ++g) {
            os << "[" << slice0(0, g) << ","
                      << slice0(1, g) << ","
                      << slice0(2, g) << "], ";
        }
        LOG(DEBUG) << "g-mesh for 1st k-point:\n" << os.str();
#endif
    }
    H5Tclose(dtype);
    H5Dclose(dset);
}
/** \brief Read the potential for the radial solver, 1 array (column) per atom type
 */
void FleurData::hdf_potential(const std::string loc)
{
    potential.set_size(max_mesh_pts, num_types);
    hid_t dset = H5Dopen(h5file, loc.c_str(), H5P_DEFAULT);
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, potential.memptr());
    H5Dclose(dset);
#ifndef NDEBUG
    std::ostringstream os;
    for (auto el : potential) {
        os << el << ", ";
    }
    LOG(DEBUG) << "Radial potential:\n" << os.str();
#endif
}

/** \brief Read the energy parameters, one per l-value & type
 *
 * "l-value" refers to the so-called spherical cutoffs, see hdf_lmax().
 */
void FleurData::hdf_energies(const std::string loc)
{
    // store as rectangular array, though some types might only have fewer valid
    // energies. needs to be handled in the getter.
    //LOG(INFO) << "Energies " << lmax;
    //LOG(INFO) << "Energies " << lmax.max()+1 << " " << num_types;
    energies.set_size(lmax.max()+1, num_types);
    hid_t dset = H5Dopen(h5file, loc.c_str(), H5P_DEFAULT);
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, energies.memptr());
    H5Dclose(dset);
    LOG(DEBUG) << "Radial energies:\n" << energies;
}
/** \brief Read per-type cutoffs of the l-parameter
 *
 * There are two different cutoffs in use in Fleur, the "spherical" and
 * "nonspherical" one. TODO: description?
 */
void FleurData::hdf_lmax(const std::string loc_l, const std::string loc_lnonsph)
{
    // one entry per atom type
    lmax.set_size(num_types);
    lmax_nonsph.set_size(num_types);
    {
        // "spherical"
        hid_t dset = H5Dopen(h5file, loc_l.c_str(), H5P_DEFAULT);
        H5Dread(dset, h5_int_mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, lmax.memptr());
        //H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, lmax.memptr());
        //H5Dread(dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, lmax.memptr()); // DIEGO
        //LOG(INFO) << "LMAX" << lmax;
        H5Dclose(dset);
    }
    {
        // "nonspherical"
        hid_t dset = H5Dopen(h5file, loc_lnonsph.c_str(), H5P_DEFAULT);
        H5Dread(dset, h5_int_mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, lmax_nonsph.memptr());
        //H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, lmax_nonsph.memptr());
        //H5Dread(dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, lmax_nonsph.memptr()); // DIEGO
        H5Dclose(dset);
    }
    LOG(DEBUG) << "l-cutoffs:\n spherical:\n" << lmax
               << "nonspherical:\n" << lmax_nonsph;
}
/** \brief Read per-type radial mesh properties
 *
 * See rad_mesh for a description of what those are.
 */
void FleurData::hdf_radmesh(const std::string locid)
{
    hid_t d_mesh = H5Dopen(h5file, locid.c_str(), H5P_DEFAULT);
    hid_t dtype = H5Dget_type(d_mesh);
    // we already know that one mesh is stored per type, so we don't
    // query the dataspace for its extent
    mesh_params.resize(num_types);
    H5Dread(d_mesh, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mesh_params[0]);
    H5Tclose(dtype);
    H5Dclose(d_mesh);
    H5LTget_attribute_int(h5file, locid.c_str(), "max_mesh_pts", &max_mesh_pts);
    std::ostringstream ss;
    for (auto mesh : mesh_params) {
        ss << mesh.num << ", " << mesh.dx << ", " << mesh.rmt << '\n';
    }
    LOG(DEBUG) << "Radial mesh:\n" << ss.str() << "max(num_pts) = " << max_mesh_pts;
}

/** \brief Read atom properties: types, positions, symmetry
 *
 * Symmetry in this case means the rotation matrices needed to rotate into the
 * representative atom for which the potential is stored.
 */
void FleurData::hdf_atoms(const std::string locid)
{
    hid_t d_atoms = H5Dopen(h5file, locid.c_str(), H5P_DEFAULT);
    H5LTget_attribute_int(h5file, "/atoms", "num_atoms", &num_atoms);
    LOG(DEBUG) << "num_atoms = " << num_atoms;
    H5LTget_attribute_int(h5file, "/atoms", "num_types", &num_types);
    LOG(DEBUG) << "num_types = " << num_types;
    types.set_size(num_atoms);
    positions.set_size(3, num_atoms);
    rotmats.set_size(3, 3, num_atoms);
    // read fields individually
    {
        // ---- atom types ----
        //hid_t t_field = H5Tcreate(H5T_COMPOUND, sizeof(int));
	hid_t t_field = H5Tcreate(H5T_COMPOUND, sizeof(arma::sword));
        //hid_t t_field = H5Tcreate(H5T_COMPOUND, sizeof(long int)); // DIEGO
        //H5Tinsert(t_field, "type", 0, H5T_NATIVE_INT);
        H5Tinsert(t_field, "type", 0, h5_int_mem_type);
        H5Dread(d_atoms, t_field, H5S_ALL, H5S_ALL, H5P_DEFAULT, types.memptr());
        //H5Dread(d_atoms, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, types.memptr()); // DIEGO
        LOG(DEBUG) << "types:\n" << types;
        //LOG(INFO) << "types:\n" << types;
        H5Tclose(t_field);
    }
    {
        // ---- positions ----
        hsize_t dims[] = {3};
        hid_t t_arr = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
        hid_t t_field = H5Tcreate(H5T_COMPOUND, 3*sizeof(double));
        H5Tinsert(t_field, "position", 0, t_arr);
        H5Dread(d_atoms, t_field, H5S_ALL, H5S_ALL, H5P_DEFAULT, positions.memptr());
        LOG(DEBUG) << "positions:\n" << positions;
        H5Tclose(t_arr);
        H5Tclose(t_field);
    }
    {
        // ---- rotation matrices ----
        hsize_t dims[] = {3,3};
        hid_t t_arr = H5Tarray_create2(H5T_NATIVE_DOUBLE, 2, dims);
        hid_t t_field = H5Tcreate(H5T_COMPOUND, 9*sizeof(double));
        H5Tinsert(t_field, "rotation", 0, t_arr);
        H5Dread(d_atoms, t_field, H5S_ALL, H5S_ALL, H5P_DEFAULT, rotmats.memptr());
        LOG(DEBUG) << "Rotation matrices (w/o rec basis):\n" << rotmats;
        H5Tclose(t_arr);
        H5Tclose(t_field);
        // apply basis matrix
        for (size_t atom = 0; atom < rotmats.n_slices; ++atom) {
            rotmats.slice(atom) *= rec_basis;
        }
    }
    H5Dclose(d_atoms);
}

/** \brief Read rotation matrices for equivalent atoms
 *
 * Fleur considers atoms of the same type that can be mapped into each other
 * through a group operation equivalent, since the potential inside the MT
 * is the same. Only when matching at the MT boundary, a rotation of the
 * coordinate system has to be applied.
 * The rotation matrix is read in the routine for each atom.
 * TODO: design goal difference, do we need equivalent atoms? probably not,
 * but this will only happen once the potential etc. is generated independently
 * of fleur.
 */
void FleurData::read_rotmat(const std::string path_prefix)
{
    RealBinmat rmat;
    rotmats.set_size(3, 3, num_atoms); // one 3x3 matrix per atom
    for (size_t atom = 1; atom <= num_atoms; ++atom) {
        std::ostringstream ss;
        ss << path_prefix << std::setw(3) << std::setfill('0') << atom << ".bin";
        try {
            rmat.read(ss.str());
            // differs from fleur: already premultiply the basis
            rotmats.slice(atom - 1) = rmat.get_data()*rec_basis;
        } catch(std::ifstream::failure e) {
            // might not be present in this data set
            LOG(ERROR) << "Could not read rotation matrix " << ss.str()
                       << ". Defaulting to Identity.";
            rotmats.slice(atom -1) = arma::eye(3,3)*rec_basis;
        }
    }
}

/** \brief Read four T-matrices calculated by Fleur in tlmplm
 *
 * The matrices read here are what Fleur considern the T-mats, not what is in
 * the formulae:
 *
 * 1. They do not contain the diagonal terms
 * 2. They contain the factor i^(l-l')
 *
 * CAUTION: Since the complex-i factor is already in, it must not appear in the
 * A-B-tensors themselves, contrary to their mathematical definition.
 */
/*
void FleurData::read_tmat(const std::string path_prefix)
{
    // fleur "terminology": u/d stands for function/derivative
    // those matrices appear in products of the A, B coefficients:
    // A^H*tuu*A, B^H*tdu*A, and so on
    std::string names[4] = {"uu", "dd", "ud", "du"},
                suffix = ".bin";
    // Read first one separately to get size
    ComplexBinmat cbm(path_prefix + names[0] + suffix);
    tmat.set_size(cbm.get_size().first, cbm.get_size().second, 4);
    tmat.slice(0) = cbm.get_data();
    assert(tmat.n_rows == tmat.n_cols); // should be square
    // read remaining matrices
    for (int i = 1; i < 4; ++i) {
        cbm.read(path_prefix + names[i] + suffix);
        tmat.slice(i) = cbm.get_data();
    }
    // each slice so far only contains the lower triangle of what are in total
    // 3 independent matrices:
    // AA and BB only need to be expanded
    tmat.slice(T_AA) = arma::symmatl(tmat.slice(T_AA));
    tmat.slice(T_BB) = arma::symmatl(tmat.slice(T_BB));
    // AB and BA need to be combined in a single matrix
    arma::cx_vec diag_AB = tmat.slice(T_AB).diag();
    tmat.slice(T_AB) =   arma::trimatl(tmat.slice(T_AB))
                    + arma::trimatu(tmat.slice(T_BA).t());
    tmat.slice(T_AB).diag() = diag_AB;
    // Sane input data: Diagonals must match
    auto diff = arma::abs(diag_AB-tmat.slice(T_BA).diag());
    bool diag_eq = arma::as_scalar(arma::all(diff < 1e-15));
    assert(diag_eq);
    tmat.slice(T_BA) = tmat.slice(T_AB).t();
}
*/
/** \brief Read kpt data output/dumped from Fleur
 *
 * Reads all k-points and associated g-point meshes.
 * \param       kpts_path     Path to the file containing the k-points in Fleur format
 * \param       gpt_prefix    appending xxx.bin with xxx in [1, num_kpts] should yield valid g-point files
 */
void FleurData::read_kpt_dump(const std::string kpts_path, const std::string gpt_prefix)
{
    using arma::span; // for submatrix view
    /* STEP 1: the 'kpts' file */
    std::ifstream kpts(kpts_path);
    int num_pts, scale;
    kpts >> num_pts >> scale;
    LOG(DEBUG) << "Number of kpts in file: " << num_pts << ", scale factor: " << scale;
    // each row of the input file has 3 k-coords (units of rec. basis vectors) and a weight
    k_points.set_size(4, num_pts);
    // TODO: find a (more) sane way to parse lines to something usable with iostream
    std::string tmp;
    kpts >> tmp; // skip empty row
    for (int j = 0; j < num_pts; ++j) {
        for (int i = 0; i < 4; ++i) {
            kpts >> k_points(i, j);
        }
    }
    // rescale all k-points (all rows except last (4th), all columns)
    k_points(span(0, 2), span::all) /= scale;

    /* STEP 2: read the g-points */
    flapw::IntBinmat gpts;
    assert(max_gpoints > 0);
    throw "broken by changes";
    /*
    g_mesh.set_size(num_pts); // one matrix of columnwise g-pts per k-point 
    for (size_t k_index = 1; k_index <= num_pts; ++k_index) {
        std::ostringstream ss;
        ss << gpt_prefix << std::setw(3) << std::setfill('0') << k_index << ".bin";
        gpts.read(ss.str());
        g_mesh(k_index - 1) = gpts.get_data();
    }
    */
}
/** \brief Read basis vectors of real and reciprocal basis
 */
void FleurData::hdf_basis(const std::string amat, const std::string bmat)
{
    hid_t dset, dtype;
    dset = H5Dopen2(h5file, amat.c_str(), H5P_DEFAULT);
    dtype = H5Dget_type(dset);
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, real_basis.memptr());
    H5Dclose(dset);
    //real_basis = rbm.get_data();
    LOG(DEBUG) << "Real basis matrix:\n" << real_basis;
    dset = H5Dopen2(h5file, bmat.c_str(), H5P_DEFAULT);
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rec_basis.memptr());
    H5Dclose(dset);
    //rbm.read(bmat);
    //rec_basis = rbm.get_data();
    LOG(DEBUG) << "Reciprocal basis matrix:\n" << rec_basis;
}

/** \brief Read potential from files [prefix][l].bin
 */
void FleurData::read_pot_files(const std::string prefix)
{
    assert(max_mesh_pts > 0); // guaranteed by read_scalars
    potential.set_size(max_mesh_pts, lmax(0)+1);
    RealBinmat vr;
    for (int l = 0; l <= lmax(0); ++l) {
        vr.read(prefix + std::to_string(l) + ".bin");
        potential.col(l) = vr.get_data();
    }
    assert( (potential.n_rows == max_mesh_pts) && (potential.n_cols == lmax(0) + 1) );
}
/** \brief Read atom positions from the given filename
 */
void FleurData::read_atom_pos(const std::string p)
{
    LOG(DEBUG) << "Reading atom positions from " << p;
    RealBinmat pos(p);
    // TODO: copies, copies everywhere!
    positions = pos.get_data();
    LOG(DEBUG) << "Read atom positions:\n" << positions;
    assert(num_atoms > 0); // should be guaranteed by read_scalars()
    assert(positions.n_cols == num_atoms);
}

/** \brief Read scalar values and dimensions of the dump files
 *
 * Should be called first, since it establishes the sizes of other data objects
 */
void FleurData::read_scalars(const std::string p)
{
    LOG(DEBUG) << "Reading scalars from " << p;
    std::ifstream f(p);
    // size of data to follow
    f >> lmax(0) >> num_atoms >> max_gpoints;
    if ((lmax(0) < 0) || (num_atoms < 1) || (max_gpoints < 1))
        throw std::invalid_argument("lmax, max_gpoints or num_atoms out of range");
    LOG(DEBUG) << "lmax = " << lmax << ", num_atoms = " << num_atoms
               << "max_gpoints = " << max_gpoints;
    // read mesh data (1 entry per atom)
    mesh_params.reserve(num_atoms);
    max_mesh_pts = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        rad_mesh r;
        // 3 values per line
        f >> r.num >> r.dx >> r.rmt;
        mesh_params.push_back(r);
        max_mesh_pts = max_mesh_pts > r.num ? max_mesh_pts : r.num;
        LOG(DEBUG) << "atoms #" << i << ": #pts=" << r.num
                   << ", dx=" << r.dx << ", r_mt=" << r.rmt;
    }
    // read lmax+1 energies
    energies.set_size(lmax(0) + 1);
    for (size_t i = 0; i <= lmax(0); ++i) {
        f >> energies(i);
    }
    LOG(DEBUG) << "Energies:\n" << energies;
}


} // namespace flapw
