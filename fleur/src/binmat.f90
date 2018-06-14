! Write a real or complex matrix into a binary file
!
! File format: small header with 3 entries:
! 1 byte: 'C' or 'R' for real or complex
! 4 bytes: integer, # of rows
! 4 bytes: integer, # of columns
! afterwards, row*col entries in column-major format, per entry:
! if real:
!   8 bytes: double, entry value
! if complex:
!   16 bytes: two doubles, first real, then imaginary part
! To accomplish this with Fortran IO, it needs `access=stream`, which is 
! available since F2013
module binmat
use hdf5
use h5lt
implicit none
private
public :: write_mat, read_mat, write_cvec, write_tlm

integer, parameter :: dp = selected_real_kind(12)

interface write_mat
    module procedure write_mat_R, write_mat_C, write_vec_C, write_vec_R, &
                     write_3vecs_I
end interface
! produce a complex vector out of different input arguments
interface write_cvec
    module procedure write_2rvecs_as_C
end interface
interface read_mat
    module procedure read_mat_R, read_mat_C
end interface

!-------------------------
!- HDF EXPORTS
!-------------------------
public :: binmat_start, binmat_end, write_hdf_mat, write_hdf_atoms, &
          write_hdf_mesh, write_hdf_tmat, &
          prepare_hdf_kpts, add_hdf_kpt, finish_hdf_kpts, write_hdf_hs
! add a named matrix as a dataset to the top-level group
interface write_hdf_mat
    module procedure write_hdf_mat_2, write_hdf_arr_int
end interface

type hdf_data_collection
! Stores run-time information about h5 handles, etc
    integer(hid_t) :: file
    !-------
    ! IO types
    integer(hid_t) :: t_atom, t_atom_type, t_atom_pos, t_atom_rot
    integer(hid_t) :: t_mesh, t_mesh_pts, t_mesh_dx, t_mesh_rmt
    integer(hid_t) :: t_kpt, t_kpt_numg, t_kpt_k, t_kpt_gmesh
    integer(hid_t) :: cx ! complex data, two doubles
    !-------
    integer(hid_t) :: g_tmat, d_kmeshes, g_hmat, g_smat
end type

character(len=*), parameter :: hdf_name = "fleur_dump.h5"
type(hdf_data_collection) :: h5dat
contains !!!!!!!!!!!!!!!!!!!!!!!
!>>>>>>> HDF5-specific stuff <<<<<<<
subroutine binmat_start()
    integer :: err
    ! Initialize hdf5
    call h5open_f(err)
    ! open dump file
    call h5fcreate_f(hdf_name, H5F_ACC_TRUNC_F, h5dat%file, err)
    ! create special groups & datatypes
    call prepare_atoms()
    call prepare_mesh()
    call prepare_tlm()
    ! complex data type, no default in HDF
    call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/2_hsize_t/), h5dat%cx, err)
    ! group for h and s matrices
    call h5gcreate_f(h5dat%file, "/matrix_hamilton", h5dat%g_hmat, err)
    call h5gcreate_f(h5dat%file, "/matrix_overlap", h5dat%g_smat, err)
contains
    subroutine prepare_atoms()
        type fleur_atom
        ! only used to create the file layout for dumped atom data
            integer :: type_id
            real(kind=dp) :: pos(3)
            integer :: rot(3,3) ! rotate into equivalent atom
        end type
        type(fleur_atom) :: da ! dummy atom, for offsets
        integer(hid_t) :: t_pos, t_rot
        integer(size_t) :: offset
        call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/3_hsize_t/), t_pos, err)
        call h5tarray_create_f(H5T_NATIVE_INTEGER, 2, (/3_hsize_t, 3_hsize_t/), t_rot, err)
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(da), h5dat%t_atom, err)
        ! insert types with offset from start of compound type
        offset = h5offsetof(c_loc(da), c_loc(da%type_id))
        call h5tinsert_f(h5dat%t_atom, "type", &
                            offset, &
                            H5T_NATIVE_INTEGER, err)
        offset = h5offsetof(c_loc(da), c_loc(da%pos(1)))
        call h5tinsert_f(h5dat%t_atom, "position", &
                            offset, &
                            t_pos, err)
        offset = h5offsetof(c_loc(da), c_loc(da%rot(1,1)))
        call h5tinsert_f(h5dat%t_atom, "rotation", &
                            offset, &
                            t_rot, err)
        ! create compound types for memory layouts to allow separate writing
        ! of type, pos and rot. matching is done by the name each component has.
        offset = 0
        ! ---> type
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(da%type_id), h5dat%t_atom_type, err)
        call h5tinsert_f(h5dat%t_atom_type, "type", offset, H5T_NATIVE_INTEGER, err)
        ! ---> position
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(da%pos), h5dat%t_atom_pos, err)
        call h5tinsert_f(h5dat%t_atom_pos, "position", offset, t_pos, err)
        ! ---> rotation
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(da%rot), h5dat%t_atom_rot, err)
        call h5tinsert_f(h5dat%t_atom_rot, "rotation", offset, t_rot, err)
        ! created all compound type, release "local" individuals (?)
        call h5tclose_f(t_pos, err)
        call h5tclose_f(t_rot, err)
        ! save the atom type for easier reading
        call h5tcommit_f(h5dat%file, "type_atoms", h5dat%t_atom, err)
    end subroutine
    subroutine prepare_mesh()
        type mesh_param
            integer :: num_pts
            real(kind=dp) :: dx, rmt
        end type
        type(mesh_param) :: dmp ! memory layout
        integer(size_t) :: offset
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(dmp), h5dat%t_mesh, err)
        offset = h5offsetof(c_loc(dmp), c_loc(dmp%num_pts))
        call h5tinsert_f(h5dat%t_mesh, "num_points", &
                           offset, &
                           H5T_NATIVE_INTEGER, err)
        offset = h5offsetof(c_loc(dmp), c_loc(dmp%dx))
        call h5tinsert_f(h5dat%t_mesh, "dx", &
                           offset, &
                           H5T_NATIVE_DOUBLE, err)
        offset = h5offsetof(c_loc(dmp), c_loc(dmp%rmt))
        call h5tinsert_f(h5dat%t_mesh, "rmt", &
                           offset, &
                           H5T_NATIVE_DOUBLE, err)
        ! as with the atom data, create compound types for memory layout
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(dmp%num_pts), h5dat%t_mesh_pts, err)
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(dmp%dx), h5dat%t_mesh_dx, err)
        call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(dmp%rmt), h5dat%t_mesh_rmt, err)
        offset = 0
        call h5tinsert_f(h5dat%t_mesh_pts, "num_points", offset, H5T_NATIVE_INTEGER, err)
        call h5tinsert_f(h5dat%t_mesh_dx, "dx", offset, H5T_NATIVE_DOUBLE, err)
        call h5tinsert_f(h5dat%t_mesh_rmt, "rmt", offset, H5T_NATIVE_DOUBLE, err)
        ! save type to file
        call h5tcommit_f(h5dat%file, "type_mesh", h5dat%t_mesh, err)
    end subroutine
    subroutine prepare_tlm()
        ! different groups of tlmplm have different sizes, thus use separate datasets
        ! and a single top-level group
        integer(hid_t) :: gid
        call h5gcreate_f(h5dat%file, "tmat", h5dat%g_tmat, err)
    end subroutine
end subroutine

subroutine write_hdf_mat_2(dset_name, data, comment)
    character(len=*), intent(in) :: dset_name
    real(kind=dp), intent(in) :: data(:, :)
    character(len=*), intent(in), optional :: comment
    ! locals
    integer :: err
    integer, parameter :: rank = 2
    integer(hsize_t) :: dims(rank)
    dims = shape(data)
    call h5ltmake_dataset_double_f(h5dat%file, dset_name, &
                rank, dims, &
                data, err)
    if (present(comment)) then
         call h5ltset_attribute_string_f(h5dat%file, dset_name, "description", comment, err)
    end if
end subroutine
subroutine write_hdf_arr_int(dset_name, data, comment)
    character(len=*), intent(in) :: dset_name
    integer, intent(in) :: data(:)
    character(len=*), intent(in), optional :: comment
    ! locals
    integer :: err
    integer, parameter :: rank = 1
    integer(hsize_t) :: dims(rank)
    dims = shape(data)
    call h5ltmake_dataset_int_f(h5dat%file, dset_name, &
                rank, dims, &
                data, err)
    if (present(comment)) then
         call h5ltset_attribute_string_f(h5dat%file, dset_name, "description", comment, err)
    end if
end subroutine
subroutine write_hdf_hs(h, s, k_index, num)
    !character(len=*), parameter :: dset_name
    complex(kind=dp), intent(in) :: h(:), s(:)
    integer, intent(in) :: k_index, num
    !character(len=*), intent(in), optional :: comment
    ! locals
    integer :: err
    integer, parameter :: rank = 1
    integer(hsize_t) :: dims(rank)
    integer(hid_t) :: dsp, dset
    dims = num*(num+1)/2 !shape(h)
    call h5screate_simple_f(rank, dims, dsp, err)
    ! h-matrix
    print*, 'write_hdf_hs: Writing, size ', dims
    call h5dcreate_f(h5dat%g_hmat, str(k_index), h5dat%cx, dsp, dset, err)
    call h5dwrite_f(dset, h5dat%cx, c_loc(h), err)
    call h5dclose_f(dset, err)
    ! s-matrix
    call h5dcreate_f(h5dat%g_smat, str(k_index), h5dat%cx, dsp, dset, err)
    call h5dwrite_f(dset, h5dat%cx, c_loc(s), err)
    call h5dclose_f(dset, err)
    call h5sclose_f(dsp, err)
!     if (present(comment)) then
!          call h5ltset_attribute_string_f(h5dat%file, dset_name, "description", comment, err)
!     end if
end subroutine

subroutine prepare_hdf_kpts(num_kpts, max_num_gpts)
! Create datatype to store all kpts, based on maximum number of g-pts
    integer, intent(in) :: num_kpts, max_num_gpts
    ! locals
    type k_mesh_dummy
        integer :: num_gpts
        real(kind=dp) :: k_point(3)
        ! need to add size by hand:
        integer :: g_mesh ! (3, max_num_gpts)
    end type
    type(k_mesh_dummy) :: dk ! dummy for type creation
    integer(hid_t) :: t_k, t_gmesh
    integer(size_t) :: offset, sz_tmp, sz_g_mesh
    integer :: err
    call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/3_hsize_t/), t_k, err)
    call h5tarray_create_f(H5T_NATIVE_INTEGER, 2, &
            (/3_hsize_t, int(max_num_gpts, hsize_t)/), t_gmesh, err)
    sz_g_mesh = 3*max_num_gpts*c_sizeof(dk%g_mesh)
    sz_tmp = c_sizeof(dk) + sz_g_mesh
    call h5tcreate_f(H5T_COMPOUND_F, sz_tmp, h5dat%t_kpt, err)
    ! insert types into compound
    offset = h5offsetof(c_loc(dk), c_loc(dk%num_gpts))
    call h5tinsert_f(h5dat%t_kpt, "num_gpts", &
                        offset, &
                        H5T_NATIVE_INTEGER, err)
    offset = h5offsetof(c_loc(dk), c_loc(dk%k_point(1)))
    call h5tinsert_f(h5dat%t_kpt, "k_point", &
                        offset, &
                        t_k, err)
    offset = h5offsetof(c_loc(dk), c_loc(dk%g_mesh))
    call h5tinsert_f(h5dat%t_kpt, "g_mesh", &
                        offset, &
                        t_gmesh, err)
    ! create compound types for memory layout to insert the data field-wise
    call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(dk%num_gpts), h5dat%t_kpt_numg, err)
    call h5tcreate_f(H5T_COMPOUND_F, c_sizeof(dk%k_point), h5dat%t_kpt_k, err)
    call h5tcreate_f(H5T_COMPOUND_F, sz_g_mesh, h5dat%t_kpt_gmesh, err)
    offset = 0
    call h5tinsert_f(h5dat%t_kpt_numg, "num_gpts", offset, H5T_NATIVE_INTEGER, err)
    call h5tinsert_f(h5dat%t_kpt_k, "k_point", offset, t_k, err)
    call h5tinsert_f(h5dat%t_kpt_gmesh, "g_mesh", offset, t_gmesh, err)
    ! write datatype to file
    call h5tcommit_f(h5dat%file, "type_kmeshes", h5dat%t_kpt, err)
    ! created all compound type, release "local" individuals (?)
    call h5tclose_f(t_k, err)
    call h5tclose_f(t_gmesh, err)
    !---------
    ! Create the dataspace and -set
    call create_dset()
    ! Create a comment describing the data
    call h5ltset_attribute_string_f(h5dat%file, "/k_meshes", "Dataformat", &
      "Each entry contains only 3*num_gpts valid entries. "//&
      "Must skip the bogus data at the end.", &
      err)
    call h5ltset_attribute_int_f(h5dat%file, "/k_meshes", "num_kpts", &
                        (/num_kpts/), 1, err)
    call h5ltset_attribute_int_f(h5dat%file, "/k_meshes", "max_gpts", &
                        (/max_num_gpts/), 1, err)
contains
    subroutine create_dset
        integer(hsize_t) :: dims(1)
        integer(hid_t) :: dsp_kpt
        dims = num_kpts
        call h5screate_simple_f(1, dims, dsp_kpt, err)
        call h5dcreate_f(h5dat%file, "/k_meshes", h5dat%t_kpt, dsp_kpt, &
                            h5dat%d_kmeshes, err)
        call h5sclose_f(dsp_kpt, err)
    end subroutine
end subroutine

subroutine add_hdf_kpt(nk, kpt, num_gpts, g1, g2, g3)
! Add a single kpt and the associated g-mesh to the datafile
    real(kind=dp), intent(in) :: kpt(3)
    integer, intent(in) :: nk, num_gpts
    integer, intent(in), dimension(:) :: g1, g2, g3
    ! locals
    integer(hid_t) :: dsp_mem, dsp_file
    integer, dimension(3, size(g1)) :: g_tmp
    integer(hsize_t) :: offset(1)
    integer :: err
    g_tmp(1, :) = g1
    g_tmp(2, :) = g2
    g_tmp(3, :) = g3
    ! TODO:
    ! * better way than to copy g?
    ! * need dataspace fiddling similar to tmat to write only 1 kpt at a time
    ! (currently breaks for more than one)
    ! select correct offset in file dataspace
    call h5dget_space_f(h5dat%d_kmeshes, dsp_file, err)
    offset = nk
    call h5sselect_elements_f(dsp_file, H5S_SELECT_SET_F, 1, 1, &
                                offset, err)
    call h5screate_simple_f(1, (/1_hsize_t/), dsp_mem, err)
    call h5dwrite_f(h5dat%d_kmeshes, h5dat%t_kpt_numg, c_loc(num_gpts), err, dsp_mem, dsp_file)
    call h5dwrite_f(h5dat%d_kmeshes, h5dat%t_kpt_k, c_loc(kpt), err, dsp_mem, dsp_file)
    call h5dwrite_f(h5dat%d_kmeshes, h5dat%t_kpt_gmesh, c_loc(g_tmp), err, dsp_mem, dsp_file)
    call h5sclose_f(dsp_mem, err)
    call h5sclose_f(dsp_file, err)
end subroutine

subroutine finish_hdf_kpts()
! cleanup stuff
    integer :: err
    call h5dclose_f(h5dat%d_kmeshes, err)
    call h5tclose_f(h5dat%t_kpt, err)
    call h5tclose_f(h5dat%t_kpt_numg, err)
    call h5tclose_f(h5dat%t_kpt_k, err)
    call h5tclose_f(h5dat%t_kpt_gmesh, err)
end subroutine

subroutine write_hdf_atoms(num, pos, rotmats, num_types, num_eq)
! Write the per-atom data with a custom datatype that groups data belonging to a
! single atom together
    integer, intent(in) :: num ! of atoms
    real(kind=dp), intent(in) :: pos(3, num)
    integer, intent(in) :: rotmats(3,3, num)
    integer, intent(in) :: num_types, num_eq(num_types)
    ! locals
    integer(hid_t) :: dsp_atoms, dset_atoms
    integer(hsize_t) :: dims(3)
    integer :: err
    integer :: at_types(num), i, at_prev
    dims(1) = num
    ! "s"pace [extent of data, logical, i.e. rank/number]
    call h5screate_simple_f(1, dims(1), dsp_atoms, err)
    ! "d"ataset [location & physical, location within h5 file]
    call h5dcreate_f(h5dat%file, "/atoms", h5dat%t_atom, dsp_atoms, &
                        dset_atoms, err)
!    If one would actually write an array of atoms:
!     at_ptr = c_loc(atom)
!     call h5dwrite_f(dset_atoms, h5dat%t_atom, &
!                         at_ptr, err)

    ! Each property of the atoms is in a separate array, thus we don't copy the
    ! data artificially into "atom" objects, but use it directly with a different
    ! HDF "memory datatype"
    ! ATOM TYPES
    ! this one's a bit tricky: need to construct an explicit list from the implicit
    ! information how many contiguous atoms are of the same type
    at_types = -1
    at_prev = 0
    do i = 1, num_types
        at_types(at_prev+1:at_prev + num_eq(i)) = i
        at_prev = at_prev + num_eq(i)
    end do
    dims(1) = num
    call h5dwrite_f(dset_atoms, h5dat%t_atom_type, at_types, (/dims(1)/), err)
    ! POSITIONS
    dims(1:2) = shape(pos)
    call h5dwrite_f(dset_atoms, h5dat%t_atom_pos, pos, dims(1:2), err)
    ! ROTATION MATRICES
    dims = shape(rotmats)
    call h5dwrite_f(dset_atoms, h5dat%t_atom_rot, rotmats, dims, err)
    ! cleanup
    call h5dclose_f(dset_atoms, err)
    call h5sclose_f(dsp_atoms, err)
    ! METADATA
    ! store actual number of atoms and different types
    call h5ltset_attribute_int_f(h5dat%file, "/atoms", "num_atoms", (/num/), 1, err)
    call h5ltset_attribute_int_f(h5dat%file, "/atoms", "num_types", (/num_types/), 1, err)
end subroutine

subroutine write_hdf_mesh(num_types, num_points, dx, rmt)
! Write information about the radial mesh
    integer, intent(in) :: num_types
    integer, intent(in) :: num_points(num_types)
    real(kind=dp), intent(in) :: dx(num_types), rmt(num_types)
    ! locals
    integer(hid_t) :: dsp_mesh, dset_mesh
    integer(hsize_t) :: dims(1)
    integer :: err
    dims = num_types
    call h5screate_simple_f(1, dims, dsp_mesh, err)
    call h5dcreate_f(h5dat%file, "/mesh", h5dat%t_mesh, dsp_mesh, &
                        dset_mesh, err)
    ! "easy" data: everything is just a 1D array with 1 entry per dspace item
    call h5dwrite_f(dset_mesh, h5dat%t_mesh_pts, num_points, dims, err)
    call h5dwrite_f(dset_mesh, h5dat%t_mesh_dx, dx, dims, err)
    call h5dwrite_f(dset_mesh, h5dat%t_mesh_rmt, rmt, dims, err)
    ! convenience: maximum mesh point number over all types
    call h5ltset_attribute_int_f(h5dat%file, "/mesh", "max_mesh_pts", &
                                 (/ maxval(num_points) /), 1, err)
    ! cleanup
    call h5dclose_f(dset_mesh, err)
    call h5sclose_f(dsp_mesh, err)
end subroutine


subroutine write_hdf_tmat(tuu, tdd, tud, tdu, ind, lmax)
! write all T-matrices at once, one dataset per atom type
! see write_tlm() subroutine for details of the fleur format
! TODO: more than 1 spin, will probably overwrite earlier spins
    complex(kind=dp), intent(in), dimension(0:, :, :) :: tuu, tdd, tud, tdu
    integer, intent(in) :: ind(0:,0:,:,:), & !L', L, atom type, spin
                           lmax(:) ! max l-value per atom type
    ! locals
    integer :: err, num_types, num_spins, lm_max, atype
    complex(kind=dp), allocatable :: T_out(:,:) ! reconstructed matrix
    integer,parameter :: spin = 1
    integer(hsize_t) :: dims(2)
    integer(hid_t) :: tmat_type, dsp_file, dset, dsp_mem
    character(len=1), parameter :: nl = char(10)
    num_types = size(tuu,2)
    num_spins = size(tuu,3)
    if (num_spins > 1) stop('not yet supported: > 1 spin')
    !call prepare_dataset()
    do atype = 1, size(lmax)
        lm_max = (lmax(atype)+1)**2
        dims = lm_max
        ! square matrix lm_max x lm_max of type complex
        call h5tarray_create_f(h5dat%cx, 2, dims, tmat_type, err)
        ! the write target (dataset in file on disk) is of size 4 matrices
        call h5screate_simple_f(1, (/4_hsize_t/), dsp_file, err)
        call h5dcreate_f(h5dat%g_tmat, str(atype), tmat_type, dsp_file, dset, err)
        ! ...but we only transfer one matrix at a time from memory
        call h5screate_simple_f(1, (/1_hsize_t/), dsp_mem, err)
        ! write the 4 matrices after one another
        call hdf_add_single_tmat(tuu, 1)
        call hdf_add_single_tmat(tdd, 2)
        call hdf_add_single_tmat(tud, 3)
        call hdf_add_single_tmat(tdu, 4)
        ! add matrix size as attribute for easier reading
        call h5ltset_attribute_int_f(dset, ".", "n", &
                                 (/ lm_max /), 1, err)
        ! ---> cleanup
        call h5dclose_f(dset, err)
        call h5sclose_f(dsp_file, err)
        call h5sclose_f(dsp_mem, err)
        ! next atom type might have different sized t-matrices
        call h5tclose_f(tmat_type, err)
    end do
    ! add an annotation so that it's clear which of the 4 matrices is which
    call h5ltset_attribute_string_f(h5dat%file, "/tmat", "Mapping_Info", &
        "Each dataset's number directly corresponds to one atom type. "&
      //"The only data is an array with 4 elements, containing the complex, "&
      //"lower triangular T-matrices stored in full matrix form."//nl&
      //"The indices (starting from 0) correspond to the"&
      //"following objects in Fleur:"// nl &
      //"0 -> tuu"//nl//"1 -> tdd"//nl//"2 -> tud"//nl//"3 -> tdu"//nl&
      //"(NB: the Fleur objects have been expanded using the index array)", err)
contains
    subroutine hdf_add_single_tmat(t_xx, offset)
    ! really, this is just to reduce writing the lines 4 times.
    ! MUST NOT be called at any other place than where it is right now, relies
    ! on the correct state of the dataspaces dps_mem, dsp_file, and anything else really.
        complex(kind=dp), intent(in) :: t_xx(0:, :, :)
        integer :: offset
        call tlm_pack_to_full(t_xx(:, atype, spin), ind(:, :, atype, spin), &
                                lmax(atype), T_out)
        call h5sselect_elements_f(dsp_file, H5S_SELECT_SET_F, 1, 1, &
                                    (/ int(offset, hsize_t) /), err)
        call h5dwrite_f(dset, tmat_type, c_loc(T_out), err, dsp_mem, dsp_file)
    end subroutine
end subroutine

subroutine binmat_end()
    integer :: err
    ! write all remanaing data
    call h5fflush_f(h5dat%file, H5F_SCOPE_GLOBAL_F, err)
    ! close individual groups & datatypes
    ! ---> atoms
    call h5tclose_f(h5dat%t_atom, err)
    call h5tclose_f(h5dat%t_atom_type, err)
    call h5tclose_f(h5dat%t_atom_pos, err)
    call h5tclose_f(h5dat%t_atom_rot, err)
    ! ---> mesh
    call h5tclose_f(h5dat%t_mesh, err)
    call h5tclose_f(h5dat%t_mesh_pts, err)
    call h5tclose_f(h5dat%t_mesh_dx, err)
    call h5tclose_f(h5dat%t_mesh_rmt, err)
    ! ---> t-matrices
    call h5gclose_f(h5dat%g_tmat, err)
    ! ---> h,s matrices
    call h5gclose_f(h5dat%g_hmat, err)
    call h5gclose_f(h5dat%g_smat, err)
    ! ---> basic data types
    call h5tclose_f(h5dat%cx, err)
    ! close dump file
    call h5fclose_f(h5dat%file, err)
    ! close hdf interface
    call h5close_f(err)
end subroutine
!>>>>>>>> HDF5 end <<<<<<<<<<
character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
subroutine tlm_pack_to_full(T, ind, lmax, T_out)
! Transform the flat array T into a 2D matrix T_out, using information stored
! in ind
! CAUTION: Need to deallocate T_out in calling context.
    complex(kind=dp), intent(in) :: T(0:) !index
    integer, intent(in) :: ind(0:,0:), & !L', L
                           lmax ! max l-value
    complex(kind=dp), allocatable, intent(out) :: T_out(:,:) ! reconstructed matrix
    ! locals
    integer :: L, Lp, i
    integer :: lm_max ! number of entries (l,m) up to lmax
    lm_max  = (lmax+1)**2
    if (allocated(T_out)) deallocate(T_out)
    allocate(T_out(0:lm_max-1, 0:lm_max-1))
    T_out = 0;
    do L = 0, lm_max - 1
        do Lp = 0, lm_max - 1
            ! get index out of mappin table
            ! encoding: -9999 -> values is 0 (gaunt-coeff rules)
            !           i >=0 -> lower triangle
            !           i < 0 -> upper triangle
            i = ind(Lp, L)
            !print '(3I8)', L, Lp, i
            if (i /= -9999) then
                if (i >= 0) then
                    ! TODO: order of indices, which one is rows/cols
                    ! (getting it wrong means hermitian transpose)
                    T_out(Lp, L) = T(i)
                else ! < 0
                    ! only lower triangles are needed:
                    ! * t_uu is hermitian, as is t_dd
                    ! * t_du^H = t_ud and vice versa
                    ! -> only output lower triangle, fix later in c-code
                    T_out(Lp, L) = 0 !conjg(T(-i, at, spin))
                endif
            end if
        end do
    end do
end subroutine
subroutine write_tlm(T, ind, lmax, fname)
! Write out the flat array T, which contains one of the T-matrices, as a 2D
! matrix using information stored in the index array `ind`
    complex(kind=dp), intent(in) :: T(0:,:,:) !index, atom type, spin
    integer, intent(in) :: ind(0:,0:,:,:), & !L', L, atom type, spin
                           lmax(:) ! max l-value per atom type
    character(len=*), intent(in) :: fname
    ! locals
    integer :: num_types, num_spins
    complex(kind=dp), allocatable :: T_out(:,:) ! reconstructed matrix
    integer,parameter :: at = 1, spin = 1
    num_types = size(T,2)
    num_spins = size(T,3)
    if (num_types > 1) stop('not yet supported: > 1 type of atoms')
    if (num_spins > 1) stop('not yet supported: > 1 spin')
    call tlm_pack_to_full(T(:, at, spin), ind(:, :, at, spin), lmax(at), T_out)
    ! matrix is reconstructed from 1D to 2D
    call write_mat(T_out, fname)
    deallocate(T_out)
end subroutine

subroutine write_vec_R(V, fname)
    real(kind=dp), intent(in) :: V(:)
    character(len=*), intent(in) :: fname
    integer :: length
    integer :: i
    length = size(V)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'R', length, 1
    do i = 1, length
        write(10) V(i)
    end do
    close(10)
end subroutine
subroutine write_vec_C(V, fname)
    complex(kind=dp), intent(in) :: V(:)
    character(len=*), intent(in) :: fname
    integer :: length
    integer :: i
    length = size(V)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'C', length, 1
    do i = 1, length
        write(10) V(i)
    end do
    close(10)
end subroutine
subroutine write_2rvecs_as_C(v1, v2, fname)
! ugly, version with 3 vecs writes as reals, with 2 as 1 complex?
! consistency... though this is only a dumb dump lib
    real, intent(in), dimension(:) :: v1, v2
    character(len=*), intent(in) :: fname
    integer :: r, c ! size of "concatenated" matrix
    integer :: i, j
    c = size(v1)
    if (c /= size(v2)) stop "Vectors are not the same length"
    r = 1 ! fixed, 2 reals are 1 complex entry, whole thing is a colvec
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'C', c, 1
    do j = 1, c
        write(10) v1(j)
        write(10) v2(j)
    end do
    close(10)
end subroutine
subroutine write_mat_R(M, fname)
    real(kind=dp), intent(in) :: M(:, :)
    character(len=*), intent(in) :: fname
    integer :: r, c ! size of M
    integer :: i, j
    r = size(M, 1)
    c = size(M, 2)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'R', r, c
    do j = 1, c
        do i = 1, r
            write(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine write_3vecs_I(v1, v2, v3, fname)
    integer, intent(in), dimension(:) :: v1, v2, v3
    character(len=*), intent(in) :: fname
    integer :: r, c ! size of "concatenated" matrix
    integer :: i, j
    c = size(v1)
    if (c /= size(v2) .or. c /= size(v3)) stop "Vectors are not the same length"
    r = 3 ! fixed, each vector is a row
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'I', r, c
    do j = 1, c        
        write(10) v1(j)
        write(10) v2(j)
        write(10) v3(j)        
    end do
    close(10)
end subroutine
subroutine write_mat_C(M, fname)
    complex(kind=dp), intent(in) :: M(:, :)
    character(len=*), intent(in) :: fname
    integer :: r, c ! size of M
    integer :: i, j
    r = size(M, 1)
    c = size(M, 2)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'C', r, c
    do j = 1, c
        do i = 1, r
            write(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine read_mat_R(M, fname)
    real(kind=dp), intent(out), allocatable :: M(:,:)
    character(len=*), intent(in) :: fname
    character(1) :: dtype
    integer :: r, c, i, j
    open(unit=10, file=fname, form='unformatted', access='stream')
    read(10) dtype, r, c
    !print*, 'Read: ', dtype, r, c
    if (dtype /= 'R') then
        print*, 'not a real matrix, type = ', dtype
        close(10)
        stop
    endif
    allocate(M(r, c))
    do j = 1, c
        do i = 1, r
            read(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine read_mat_C(M, fname)
    complex(kind=dp), intent(out), allocatable :: M(:,:)
    character(len=*), intent(in) :: fname
    character(1) :: dtype
    integer :: r, c, i, j
    open(unit=10, file=fname, form='unformatted', access='stream')
    read(10) dtype, r, c
    !print*, 'Read: ', dtype, r, c
    if (dtype /= 'C') then
        print*, 'not a complex matrix, type = ', dtype
        close(10)
        stop
    endif
    allocate(M(r, c))
    do j = 1, c
        do i = 1, r
            read(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
end module
