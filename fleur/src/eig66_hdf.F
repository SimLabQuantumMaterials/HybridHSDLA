      MODULE m_eig66_hdf
c*****************************************************************
c DESC:Module for hdf-io of eig-file
c      To be compatible with f90 interface of HDF, use kind for vars
c
c      !ATTENTION before calling openeig and after calling closeeig!
c      !the hdf library has to be initialized or finalized, respectively !
c
c      CONTAINS the following subroutines:
c      openeig        opens file
c      closeeig       closes file
c      read_keb       reads kpt, enpara and basis data
c      read_neig      read no of eigenvalues (and eigenvalues itself)
c      read_eig       reads eigenvectors
c      writeeig       saves all data for kpt
c      writesingleeig saves data for one kpt and energy
c
c
c                          Daniel Wortmann, Tue Nov  5 12:07:52 2002
c*****************************************************************

      USE hdf5
      USE m_hdf_tools
      IMPLICIT NONE

      PRIVATE
      INTEGER, PARAMETER :: one=1,two=2,three=3,zero=0
                               !to have the correct
                               !type for array constructors
      INTEGER(HID_T), SAVE :: fid
      INTEGER(HID_T), SAVE :: elsetid,evacsetid,ellosetid
      INTEGER(HID_T), SAVE :: bksetid,wksetid,ksetid
      INTEGER(HID_T), SAVE :: neigsetid,nvsetid,nmatsetid
      INTEGER(HID_T), SAVE :: energysetid,evsetid
      PUBLIC openeig,closeeig
      PUBLIC read_eig,read_kptenpar,read_vec,read_basis
      PUBLIC read_enpar,read_vecc,read_kpt
      PUBLIC writeeig,writesingleeig,writeeigc,writebas
C 
      CONTAINS 
!----------------------------------------------------------------------
      SUBROUTINE openeig(
     >                   nmat,neig,nkpts,jspins,lmax,nlo,ntype,
     >                   create,readonly)

c*****************************************************************
c     opens hdf-file for eigenvectors+values
c*****************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nmat,neig,nkpts,jspins,nlo,ntype,lmax
      LOGICAL, INTENT(IN) :: create,readonly
      
      INTEGER         :: hdferr,access_mode
      INTEGER(HID_T)  :: creation_prp,access_prp,spaceid
      LOGICAL         :: l_exist
      INTEGER(HSIZE_T):: dims(5)

      !Set creation and access properties
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      IF (readonly) THEN
      access_prp=H5P_DEFAULT_f
      creation_prp=H5P_DEFAULT_f
      ELSE
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)
      CALL h5pset_fapl_mpio_f(access_prp, MPI_COMM_WORLD, MPI_INFO_NULL,
     +     hdferr)
      creation_prp=H5P_DEFAULT_f !no special creation property
      ENDIF
#else
      access_prp=H5P_DEFAULT_f
      creation_prp=H5P_DEFAULT_f
#endif 
      !set access_flags according 
      IF (readonly) THEN
         access_mode=H5F_ACC_RDONLY_F
      ELSE 
         access_mode=H5F_ACC_RDWR_F 
      ENDIF
c     OPEN FILE and get FID's
      IF (create) THEN
         INQUIRE(FILE='eig.hdf',EXIST=l_exist)
         access_mode=H5F_ACC_TRUNC_F
!         IF (l_exist) WRITE (*,*)'Warning: eig.hdf was overwritten'
         CALL h5fcreate_f('eig.hdf',access_Mode, fid, hdferr
     +        ,creation_prp,access_prp)
! create dataspaces and datasets
!   scalars
         dims(:2)=(/nkpts,jspins/)
         CALL h5screate_simple_f(2,dims(:2),spaceid,hdferr) 
         CALL h5dcreate_f(fid, "neig", H5T_NATIVE_INTEGER, spaceid, 
     +                 neigsetid, hdferr)
         CALL h5dcreate_f(fid, "wk", H5T_NATIVE_DOUBLE,spaceid, 
     +                 wksetid, hdferr)
         CALL h5dcreate_f(fid, "nv", H5T_NATIVE_INTEGER, spaceid, 
     +                 nvsetid, hdferr)
         CALL h5dcreate_f(fid, "nmat", H5T_NATIVE_INTEGER, spaceid, 
     +                 nmatsetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)
!   vectors
         dims(1:3)=(/two,nkpts,jspins/)
         CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr) 
         CALL h5dcreate_f(fid, "evac", H5T_NATIVE_DOUBLE, spaceid, 
     +                 evacsetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)
         dims(:3)=(/three,nkpts,jspins/)
         CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr) 
         CALL h5dcreate_f(fid, "bk", H5T_NATIVE_DOUBLE, spaceid, 
     +                 bksetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)
#ifdef CPP_SOC 
         dims(:3)=(/2*neig,nkpts,jspins/)
#else
         dims(:3)=(/neig,nkpts,jspins/)
#endif
         CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr) 
!     ew
         CALL h5dcreate_f(fid, "energy", H5T_NATIVE_DOUBLE, spaceid, 
     +                 energysetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)
!     enparas
         dims(1:4)=(/lmax+1,ntype,nkpts,jspins/)
         CALL h5screate_simple_f(4,dims(1:4),spaceid,hdferr) 
         CALL h5dcreate_f(fid, "el", H5T_NATIVE_DOUBLE, spaceid, 
     +                 elsetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)
         
         dims(:4)=(/nlo,ntype,nkpts,jspins/)
         CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr) 
         CALL h5dcreate_f(fid, "ello", H5T_NATIVE_DOUBLE, spaceid, 
     +                 ellosetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)         
!     ev 
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
         dims(:5)=(/one,nmat,neig,nkpts,jspins/)
#else
#ifdef CPP_SOC
         dims(:5)=(/two,nmat,2*neig,nkpts,jspins/)
#else
         dims(:5)=(/two,nmat,neig,nkpts,jspins/)
#endif
#endif
         CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr) 
         CALL h5dcreate_f(fid, "ev", H5T_NATIVE_DOUBLE, spaceid, 
     +                 evsetid, hdferr)
         CALL h5sclose_f(spaceid,hdferr)
!      basis
         dims(:4)=(/nmat,three,nkpts,jspins/)
         CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)          
         CALL h5dcreate_f(fid, "k", H5T_NATIVE_INTEGER, spaceid, 
     +                 ksetid, hdferr)         
         CALL h5sclose_f(spaceid,hdferr)         
      ELSE
         CALL h5fopen_f ('eig.hdf', access_Mode, fid, hdferr,access_prp)       
         !get dataset-ids
         CALL h5dopen_f(fid, 'el', elsetid, hdferr)         
         CALL h5dopen_f(fid, 'evac', evacsetid, hdferr)         
         CALL h5dopen_f(fid, 'ello', ellosetid, hdferr)         
         CALL h5dopen_f(fid, 'bk', bksetid, hdferr)
         CALL h5dopen_f(fid, 'wk', wksetid, hdferr)
         CALL h5dopen_f(fid, 'energy', energysetid, hdferr)
         CALL h5dopen_f(fid, 'k', ksetid, hdferr)
         CALL h5dopen_f(fid, 'neig', neigsetid, hdferr)
         CALL h5dopen_f(fid, 'ev', evsetid, hdferr)
         CALL h5dopen_f(fid, 'nv', nvsetid, hdferr)
         CALL h5dopen_f(fid, 'nmat', nmatsetid, hdferr)         
      ENDIF
      IF (.not.access_prp==H5P_DEFAULT_f) CALL H5Pclose_f(access_prp
     +     ,hdferr)
      
      END SUBROUTINE openeig
!----------------------------------------------------------------------
      SUBROUTINE closeeig()
c*****************************************************************
c     closes hdf-file for eigenvectors+values
c*****************************************************************      IMPLICIT NONE
      !close datasets
      INTEGER::hdferr
      CALL h5dclose_f(elsetid,hdferr)
      CALL h5dclose_f(evacsetid,hdferr)
      CALL h5dclose_f(ellosetid,hdferr)
      CALL h5dclose_f(bksetid,hdferr)
      CALL h5dclose_f(wksetid,hdferr)
      CALL h5dclose_f(energysetid,hdferr)
      CALL h5dclose_f(ksetid,hdferr)
      CALL h5dclose_f(neigsetid,hdferr)
      CALL h5dclose_f(evsetid,hdferr)
      CALL h5dclose_f(nvsetid,hdferr)
      CALL h5dclose_f(nmatsetid,hdferr)
      !close file
      CALL h5fclose_f(fid,hdferr)
      END SUBROUTINE
!---------------------------------------------------------------------- 
      SUBROUTINE read_vec(
     >                    nk,jspin,neig,n_start,n_end,
     <                    nmat,z)

      USE m_hdf_tools
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: nk,jspin
      INTEGER, INTENT(IN)  :: neig,n_start,n_end
      INTEGER, INTENT(OUT) :: nmat
#ifdef CPP_INVERSION
      REAL,    INTENT(OUT) :: z(:,:)
#ifdef CPP_SOC
      REAL, ALLOCATABLE    :: z1(:,:,:)
#endif
#else
      COMPLEX, INTENT(OUT) :: z(:,:)
      REAL, ALLOCATABLE    :: z1(:,:,:)
#endif
      INTEGER i,j,neig_l

      neig_l = n_end - n_start + 1

! read matrix size
      CALL io_read_integer0(
     >                      nmatsetid,(/nk,jspin/),(/1,1/),
     <                                                nmat)

      IF ( nmat > size(z,1) .OR. neig_l > size(z,2) ) THEN
         WRITE (6,*) nmat,size(z,1),neig,size(z,2)
         STOP 'eig66_hdf$read_vec'
      ENDIF

      !read eigenvectors
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      CALL io_read_real2(evsetid,(/1,1,n_start,nk,jspin/),
     >                           (/1,nmat,neig_l,1,1/),
     <                           z(:nmat,:neig_l) )
#else
      ALLOCATE (z1(2,nmat,neig_l))
      CALL io_read_real3(evsetid,(/1,1,n_start,nk,jspin/),
     >                       (/2,nmat,neig_l,1,1/),z1)

      DO i=1,neig_l
         DO j=1,nmat 
#ifndef CPP_INVERSION
            z(j,i) = cmplx( z1(1,j,i) ,z1(2,j,i) )
#else
            z(j,i) = z1(1,j,i)
#endif
         ENDDO
      ENDDO

      DEALLOCATE (z1)
#endif

      END SUBROUTINE read_vec

!--------------------------------------------------------------------------
      SUBROUTINE read_kptenpar(
     >                         nk,jspin,
     <                         bk,wk,el,ello,evac)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: nk,jspin
      REAL,    INTENT(OUT) :: bk(3),evac(2),ello(:,:),el(:,:),wk

!   read enparas
      CALL io_read_real2(elsetid,(/1,1,nk,jspin/),
     +                   (/size(el,1),size(el,2),1,1/),el)
      CALL io_read_real2(ellosetid,(/1,1,nk,jspin/),
     +                   (/size(ello,1),size(ello,2),1,1/),ello)
      CALL io_read_real1(evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)

!   read kpts
      CALL io_read_real1(bksetid,(/1,nk,jspin/),(/3,1,1/),bk)
      CALL io_read_real0(wksetid,(/nk,jspin/),(/1,1/),wk)

      END SUBROUTINE read_kptenpar
!--------------------------------------------------------------------------
      SUBROUTINE read_kpt(
     >                    nk,jspin,
     <                    bk,wk)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: nk,jspin
      REAL,    INTENT(OUT) :: bk(3),wk

!   read kpts
      CALL io_read_real1(bksetid,(/1,nk,jspin/),(/3,1,1/),bk)
      CALL io_read_real0(wksetid,(/nk,jspin/),(/1,1/),wk)

      END SUBROUTINE read_kpt
!--------------------------------------------------------------------------
      SUBROUTINE read_enpar(
     >                      nk,jspin,
     <                      el,ello,evac)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: nk,jspin
      REAL,    INTENT(OUT) :: evac(2),ello(:,:),el(:,:)

!   read enparas
      CALL io_read_real2(elsetid,(/1,1,nk,jspin/),
     +                   (/size(el,1),size(el,2),1,1/),el(:,:))
      CALL io_read_real2(ellosetid,(/1,1,nk,jspin/),
     +                   (/size(ello,1),size(ello,2),1,1/),ello(:,:))
      CALL io_read_real1(evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)

      END SUBROUTINE read_enpar
!------------------------------------------------------------------------
      SUBROUTINE read_basis(nk,jspin,nlotot,
     <                      nmat,nv,k1,k2,k3,kveclo)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: nk,jspin,nlotot
      INTEGER, INTENT(OUT) :: nv,nmat
      INTEGER, INTENT(OUT) :: k1(:),k2(:),k3(:),kveclo(nlotot)

      INTEGER k,k1_t,k2_t,k3_t,n

      !read dims
      CALL io_read_integer0(nvsetid,(/nk,jspin/),(/1,1/),nv)
      CALL io_read_integer0(nmatsetid,(/nk,jspin/),(/1,1/),nmat)
      if (nv>size(k1)) STOP 'eig66_hdf$read_basis'
      !read basis
      CALL io_read_integer1(ksetid,(/1,1,nk,jspin/),(/nv,1,1,1/),k1(:nv)
     +     )
      CALL io_read_integer1(ksetid,(/1,2,nk,jspin/),(/nv,1,1,1/),k2(:nv)
     +     )
      CALL io_read_integer1(ksetid,(/1,3,nk,jspin/),(/nv,1,1,1/),k3(:nv)
     +     )
      DO k = 1, nlotot
        CALL io_read_integer0(
     >        ksetid,(/nv+k,1,nk,jspin/),(/1,1,1,1/),k1_t)
        CALL io_read_integer0(
     >        ksetid,(/nv+k,2,nk,jspin/),(/1,1,1,1/),k2_t)
        CALL io_read_integer0(
     >        ksetid,(/nv+k,3,nk,jspin/),(/1,1,1,1/),k3_t)
        DO n = 1, nv
         IF (( (k1_t == k1(n)).AND.(k2_t == k2(n)) ).AND.
     +         (k3_t == k3(n)) ) THEN
            kveclo(k) = n
            CYCLE
         ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE read_basis
!------------------------------------------------------------------------
      SUBROUTINE read_eig(
     >                    nk,jspin,
     <                    neig,eig)
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nk,jspin
      INTEGER, INTENT (OUT) :: neig
      REAL,    INTENT (OUT),OPTIONAL :: eig(:)
                                                   ! read No of eigenvalues

      CALL io_read_integer0(neigsetid,(/nk,jspin/),(/1,1/),neig)

      IF ( present(eig) ) THEN                           ! read eigenvalues
         IF ( neig > size(eig) ) THEN
           WRITE(*,*) neig,size(eig)
           STOP 'eig66_hdf$readeig'
         ENDIF
         CALL io_read_real1(energysetid,(/1,nk,jspin/),(/neig,1,1/),
     <                      eig(:neig))
      ENDIF

      END SUBROUTINE read_eig

!-------------------------------------------------------------------------

      SUBROUTINE writeeig(
     >                    nk,jspin,neig,nv,nmat,k1,k2,k3,bk,wk,
     >                    eig,el,ello,evac,z,
     >                    nlotot,kveclo,n_size,n_rank,neigd)

c*****************************************************************
c     writes all eignevecs for the nk-th kpoint
c*****************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nk,jspin,n_size,n_rank,neigd
      REAL,    INTENT(IN) :: wk
      INTEGER, INTENT(IN) :: neig,nv,nmat,nlotot
      INTEGER, INTENT(IN) :: k1(:),k2(:),k3(:),kveclo(nlotot)
      REAL,    INTENT(IN) :: bk(3),eig(:),el(:,:),evac(2),ello(:,:)
#ifdef CPP_INVERSION
      REAL,    INTENT(IN)  :: z(:,:)
#ifdef CPP_SOC
      REAL,    ALLOCATABLE :: z1(:,:,:)
#endif
#else
      COMPLEX, INTENT(IN)  :: z(:,:)
      REAL,    ALLOCATABLE :: z1(:,:,:)
#endif
      INTEGER i,j,k

!
!write enparas
!
      CALL io_write_real2(
     >                    elsetid,(/1,1,nk,jspin/),
     >                    (/size(el,1),size(el,2),1,1/),el)
      CALL io_write_real2(
     >                    ellosetid,(/1,1,nk,jspin/),
     >                    (/size(ello,1),size(ello,2),1,1/),ello)
      CALL io_write_real1(
     >                    evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)
!
!write kpts
!
      CALL io_write_real1(
     >                    bksetid,(/1,nk,jspin/),(/3,1,1/),bk)
      CALL io_write_real0(
     >                    wksetid,(/nk,jspin/),(/1,1/),wk)
!
!write basis
!
      CALL io_write_integer0(
     >                       nvsetid,(/nk,jspin/),(/1,1/),nv)
      CALL io_write_integer0(
     >                       nmatsetid,(/nk,jspin/),(/1,1/),nmat)
      CALL io_write_integer1(
     >              ksetid,(/1,1,nk,jspin/),(/nv,1,1,1/),k1(:nv))
      CALL io_write_integer1(
     >              ksetid,(/1,2,nk,jspin/),(/nv,1,1,1/),k2(:nv))
      CALL io_write_integer1(
     >              ksetid,(/1,3,nk,jspin/),(/nv,1,1,1/),k3(:nv))
      DO k = 1, nlotot
        CALL io_write_integer0(
     >        ksetid,(/nv+k,1,nk,jspin/),(/1,1,1,1/),k1(kveclo(k)))
        CALL io_write_integer0(
     >        ksetid,(/nv+k,2,nk,jspin/),(/1,1,1,1/),k2(kveclo(k)))
        CALL io_write_integer0(
     >        ksetid,(/nv+k,3,nk,jspin/),(/1,1,1,1/),k3(kveclo(k)))
      ENDDO
!
!write eigenvalues
!
      IF ( n_size > 1 ) THEN
        CALL io_write_integer0(
     >                         neigsetid,(/nk,jspin/),(/1,1/),neigd)
      ELSE
        CALL io_write_integer0(
     >                         neigsetid,(/nk,jspin/),(/1,1/),neig)
      ENDIF
      CALL io_write_real1s(
     >                     energysetid,(/n_rank+1,nk,jspin/),       ! use n_rank as offsett
     >                     (/neig,1,1/),eig(:neig),(/n_size,1,1/))  ! and n_size as stride
!                                                                     for ev-parallelization
!write eigenvectors
!
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      CALL io_write_real2s(
     >                     evsetid,(/1,1,n_rank+1,nk,jspin/),
     >            (/1,nmat,neig,1,1/),z(:nmat,:neig),(/1,1,n_size,1,1/))
#else
      ALLOCATE (z1(2,nmat,neig))

      DO i=1,neig
         DO j=1,nmat
            z1(1,j,i) = real(z(j,i))
#ifndef CPP_INVERSION
            z1(2,j,i) = aimag(z(j,i))
#else
            z1(2,j,i) = 0.0
#endif
         ENDDO
      ENDDO

!      CALL io_write_real3(evsetid,(/1,1,1,nk,jspin/),
!     >                        (/2,nmat,neig,1,1/),z1)
      CALL io_write_real3s(evsetid,(/1,1,n_rank+1,nk,jspin/),
     >                       (/2,nmat,neig,1,1/),z1,(/1,1,n_size,1,1/))
      DEALLOCATE (z1)
#endif
      END SUBROUTINE writeeig

!-------------------------------------------------------------------------

      SUBROUTINE writesingleeig(
     >                          nk,jspin,nstart,nstop,nv,nmat,eig,z,
     >                          k1,k2,k3,bk,wk,el,ello,evac)

c*****************************************************************
c     writes all eignevecs for the nk-th kpoint
c*****************************************************************
      IMPLICIT NONE
      INTEGER,INTENT(IN)         ::nk,jspin,nstart,nstop,nmat,nv
      REAL,INTENT(IN),OPTIONAL   ::eig(:)
      INTEGER,INTENT(IN),OPTIONAL::k1(:),k2(:),k3(:)
      REAL,INTENT(IN),OPTIONAL   ::bk(3),el(:,:),evac(2),
     +     ello(:,:),wk
#ifdef CPP_INVERSION
      REAL,    INTENT(IN)   ::z(:,:)
#else
      COMPLEX, INTENT(IN)   ::z(:,:)
#endif
      INTEGER ::neig

      !write enparas
      IF (present(el)) THEN
         CALL io_write_real2(elsetid,(/1,1,nk,jspin/),(/size(el,1)
     +        ,size(el,2),1,1/),el)
         CALL io_write_real2(ellosetid,(/1,1,nk,jspin/),(/size(ello,1)
     +        ,size(ello,2),1,1/),ello)
         CALL io_write_real1(evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)
      ENDIF
      IF (present(bk)) THEN !assumes that wk is also present!
      !write kpts
         CALL io_write_real1(bksetid,(/1,nk,jspin/),(/3,1,1/),bk)
         CALL io_write_real0(wksetid,(/nk,jspin/),(/1,1/),wk)
      ENDIF
      !write basis
      CALL io_write_integer0(nvsetid,(/nk,jspin/),(/1,1/),nv)
      CALL io_write_integer0(nmatsetid,(/nk,jspin/),(/1,1/),nmat)
      IF (present(k1)) THEN !assumes that all of k1..k3 are present
         CALL io_write_integer1(ksetid,(/1,1,nk,jspin/),(/nv,1,1,1/)
     +        ,k1(:nv))
         CALL io_write_integer1(ksetid,(/1,2,nk,jspin/),(/nv,1,1,1/)
     +        ,k2(:nv))
         CALL io_write_integer1(ksetid,(/1,3,nk,jspin/),(/nv,1,1,1/)
     +        ,k3(:nv))
      ENDIF
      !write eigenvalues
      IF (nstop>nstart-1) THEN ! IF nstop<nstart nothing will be written!!!
         CALL io_read_integer0(neigsetid,(/nk,jspin/),(/1,1/),neig)
         !In parallel execution this part might be wrong!
         !race condition might occur!!
         neig=max(neig,nstop-nstart+1) !determine new maximum!!
         CALL io_write_integer0(neigsetid,(/nk,jspin/),(/1,1/),neig)
         CALL io_write_real1(energysetid,(/nstart,nk,jspin/),(/nstop
     +        -nstart+1,1,1/),eig(:nstop))
         !write vectors
#ifdef CPP_INVERSION
         CALL io_write_real2(evsetid,(/1,1,nstart,nk,jspin/),(/1,nmat
     +        ,nstop-nstart+1,1,1/),z(:nmat,:nstop))
#else
         CALL io_write_complex2(evsetid,(/-1,1,nstart,nk,jspin/),(/1
     +        ,nmat,nstop-nstart+1,1,1/),z(:nmat,:nstop))
#endif
      ENDIF
      END SUBROUTINE writesingleeig

!-------------------------------------------------------------------------

      SUBROUTINE writeeigc(
     >                    nk,jspin,neig,nmat,
     >                    eig,z)

c*****************************************************************
c     writes all eignevecs for the nk-th kpoint 
c     as above but only complex
c*****************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nk,jspin
      INTEGER, INTENT(IN) :: neig,nmat
      REAL,    INTENT(IN) :: eig(:)
      COMPLEX, INTENT(IN) :: z(:,:)
      INTEGER i,j
      REAL, ALLOCATABLE :: z1(:,:,:)
!
!write eigenvalues
!
      CALL io_write_integer0(
     >                       neigsetid,(/nk,jspin/),(/1,1/),neig)
      CALL io_write_real1(
     >                    energysetid,(/1,nk,jspin/),
     >                    (/neig,1,1/),eig(:neig)) 
!                                                    
!write eigenvectors
!
      ALLOCATE (z1(2,nmat,neig))

      DO i=1,neig
         DO j=1,nmat
            z1(1,j,i) = real(z(j,i))
#ifdef CPP_INVERSION
            z1(2,j,i) = aimag(z(j,i))
#else
            z1(2,j,i) = aimag(z(j,i))
#endif
         ENDDO
      ENDDO

      CALL io_write_real3(evsetid,(/1,1,1,nk,jspin/),
     >                        (/2,nmat,neig,1,1/),z1)

      DEALLOCATE (z1)

      END SUBROUTINE writeeigc
!----------------------------------------------------------------------
      SUBROUTINE read_vecc(
     >                     nk,jspin,neig,n_start,n_end,
     <                     nmat,z)

      USE m_hdf_tools
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: nk,jspin
      INTEGER, INTENT(IN)  :: neig,n_start,n_end
      INTEGER, INTENT(OUT) :: nmat
      COMPLEX, INTENT(OUT) :: z(:,:)

      REAL, ALLOCATABLE :: z1(:,:,:)
      INTEGER i,j,neig_l

      neig_l = n_end - n_start + 1

! read matrix size
      CALL io_read_integer0(
     >                      nmatsetid,(/nk,jspin/),(/1,1/),
     <                                                nmat)

      IF ( nmat > size(z,1) .OR. neig_l > size(z,2) ) THEN
         WRITE (6,*) nmat,size(z,1),neig,size(z,2)
         STOP 'eig66_hdf$read_vec'
      ENDIF

! read eigenvectors
      ALLOCATE (z1(2,nmat,neig_l))
      CALL io_read_real3(evsetid,(/1,1,n_start,nk,jspin/),
     >                      (/2,nmat,neig_l,1,1/),z1)
      
      DO i=1,neig_l
         DO j=1,nmat
            z(j,i) = cmplx( z1(1,j,i) ,z1(2,j,i) )
         ENDDO
      ENDDO

      DEALLOCATE (z1)

      END SUBROUTINE read_vecc
!-------------------------------------------------------------------------

      SUBROUTINE writebas(
     >                    nk,jspin,nv,k1,k2,k3,nlotot,kveclo,
     >                    el,ello,evac)
c*****************************************************************
c     writes all eignevecs for the nk-th kpoint
c*****************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nk,jspin
      INTEGER, INTENT(IN) :: nv,nlotot
      INTEGER, INTENT(IN) :: k1(:),k2(:),k3(:),kveclo(nlotot)
      REAL,    INTENT(IN) :: el(:,:),evac(2),ello(:,:)

      INTEGER k
!
!write enparas
!
      CALL io_write_real2(
     >                    elsetid,(/1,1,nk,jspin/),
     >                    (/size(el,1),size(el,2),1,1/),el)
      CALL io_write_real2(
     >                    ellosetid,(/1,1,nk,jspin/),
     >                    (/size(ello,1),size(ello,2),1,1/),ello)
      CALL io_write_real1(
     >                    evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)
!
!write basis
!
      CALL io_write_integer0(
     >                       nvsetid,(/nk,jspin/),(/1,1/),nv)
      CALL io_write_integer1(
     >              ksetid,(/1,1,nk,jspin/),(/nv,1,1,1/),k1(:nv))
      CALL io_write_integer1(
     >              ksetid,(/1,2,nk,jspin/),(/nv,1,1,1/),k2(:nv))
      CALL io_write_integer1(
     >              ksetid,(/1,3,nk,jspin/),(/nv,1,1,1/),k3(:nv))
      DO k = 1, nlotot
        CALL io_write_integer0(
     >        ksetid,(/nv+k,1,nk,jspin/),(/1,1,1,1/),k1(kveclo(k)))
        CALL io_write_integer0(
     >        ksetid,(/nv+k,2,nk,jspin/),(/1,1,1,1/),k2(kveclo(k)))
        CALL io_write_integer0(
     >        ksetid,(/nv+k,3,nk,jspin/),(/1,1,1,1/),k3(kveclo(k)))
      ENDDO
!
      END SUBROUTINE writebas


      END MODULE

