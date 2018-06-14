      MODULE m_mpi_col_dos
!
! collects the information necessary for the determination of the
! DOS  from the nodes and writes it to the tmp_dos (and tmp_Ek)
! pe 0; only for parallel applications without direct global access
! files.                                                        gb 05
!     
! Now also for eigenvalue parallelization                       gb 07
!
      CONTAINS
      SUBROUTINE mpi_col_dos(
     >                       l_evp,noccbd,noccbd_l,
     >                       irank,isize,jspd,neigd,lmaxd,ntypd,
     >                       nsld,natd,ncored,nstd,nstars,layerd,
     >                       l_mcd,ndir,nkpt,jspin,ikpt,
     >                       bkpt,wk,nbands,eig,qal,qvac,qis,
     >                       qvlay,qstars,ksym,jsym,mcd,
     >                       qintsl,qmtsl,qmtp,orbcomp)

      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT (IN) :: noccbd,noccbd_l
      INTEGER, INTENT (IN) :: jspd,neigd,lmaxd,ntypd
      INTEGER, INTENT (IN) :: nsld,natd,ncored
      INTEGER, INTENT (IN) :: nstd,nstars,layerd,ndir,nkpt
      INTEGER, INTENT (IN) :: irank,isize,jspin,ikpt
      LOGICAL, INTENT (IN) :: l_mcd,l_evp
      INTEGER, INTENT (INOUT) :: nbands
      REAL,    INTENT (INOUT) :: wk

      COMPLEX, INTENT (INOUT) :: qstars(nstars,neigd,layerd,2)
      COMPLEX, INTENT (INOUT) :: mcd(3*ntypd,nstd,neigd)
      INTEGER, INTENT (INOUT) :: ksym(neigd),jsym(neigd)
      REAL,    INTENT (INOUT) :: eig(neigd),bkpt(3)
      REAL,    INTENT (INOUT) :: qis(neigd)
      REAL,    INTENT (INOUT) :: qal(0:3,ntypd,neigd)
      REAL,    INTENT (INOUT) :: qvac(neigd,2)
      REAL,    INTENT (INOUT) :: qvlay(neigd,layerd,2)
      REAL,    INTENT (INOUT) :: qintsl(nsld,neigd)
      REAL,    INTENT (INOUT) :: qmtsl(nsld,neigd)
      REAL,    INTENT (INOUT) :: orbcomp(neigd,23,natd)
      REAL,    INTENT (INOUT) :: qmtp(neigd,natd)

      INTEGER i,ie,is,id,ir,n_start,n_end,nbdim,n_rd,noccb1
      INTEGER ierr(3),stt(MPI_STATUS_SIZE)
      REAL, ALLOCATABLE :: all(:)

      IF (l_evp) THEN
       nbdim = noccbd
      ELSE
       nbdim = neigd
      ENDIF

      id = 5 + nbdim*(6 + 4*ntypd + layerd*2*(1 + nstars))
      ALLOCATE ( all(0:id) )
      IF (irank.ne.0) THEN
         all(0)       = nkpt*(jspin-1)+ikpt
         all(1:3)     = bkpt(1:3)
         all(4)       = wk
         all(5)       = real(noccbd) ! real(nbands)
         is = 6  ; ie = is - 1 + nbdim 
         all(is:ie)   = eig(1:nbdim)
         is =ie+1; ie = is - 1 + 4*ntypd*nbdim
         all(is:ie)   = RESHAPE(qal(:,:,1:nbdim),(/4*ntypd*nbdim/))
         is =ie+1; ie = is - 1 + 2*nbdim
         all(is:ie)   = RESHAPE(qvac(1:nbdim,:),(/2*nbdim/))
         is =ie+1; ie = is - 1 + nbdim
         all(is:ie)   = qis(1:nbdim)
         is =ie+1; ie = is - 1 + nbdim*layerd*2
         all(is:ie)   = RESHAPE(qvlay(1:nbdim,:,:),(/nbdim*layerd*2/))
         is =ie+1; ie = is - 1 + nstars*nbdim*layerd*2
         all(is:ie)=RESHAPE(qstars(:,1:nbdim,:,:),
     +                            (/nstars*nbdim*layerd*2/))
     
         is =ie+1; ie = is - 1 + nbdim
         all(is:ie)   = real(jsym(1:nbdim))
         is =ie+1; ie = is - 1 + nbdim
         IF (ie.ne.id) STOP 'error in mpi_col_dos!'
         all(is:ie)   = real(ksym(1:nbdim))
         CALL MPI_SEND(all,ie,MPI_DOUBLE_PRECISION,0,irank,
     +                                 MPI_COMM_WORLD,ierr)
      ELSE
         n_rd = isize-1
         IF (.not.l_evp) n_rd = min(nkpt-ikpt,isize-1)

         DO i = 1,n_rd
           CALL MPI_RECV(all,id,MPI_DOUBLE_PRECISION,i,i,
     +                           MPI_COMM_WORLD,stt,ierr)

            ir        = NINT(all(0))
            bkpt(1:3) = all(1:3)
            wk        = all(4)
            noccb1 = NINT(all(5))               ! nbands    = NINT(all(5))
            IF (l_evp) THEN                     ! collect bands
              n_start = i*noccbd_l + 1
              n_end = n_start + noccb1 - 1

              is = 6  ; ie = is - 1 + noccb1 ! bdim
              eig(n_start:n_end) = all(is:ie)
              is =ie+1; ie = is - 1 + 4*ntypd*noccb1 ! bdim
              qal(:,:,n_start:n_end) = 
     +                          RESHAPE(all(is:ie),(/4,ntypd,noccb1/))
              is =ie+1; ie = is - 1 + 2*noccb1 ! bdim
              qvac(n_start:n_end,:) = RESHAPE(all(is:ie),(/noccb1,2/))
              is =ie+1; ie = is - 1 + noccb1 ! bdim
              qis(n_start:n_end) = all(is:ie)
              is =ie+1; ie = is - 1 + noccb1*layerd*2
              qvlay(n_start:n_end,:,:) = 
     +                         RESHAPE(all(is:ie),(/noccb1,layerd,2/))
              is =ie+1; ie = is - 1 + nstars*noccb1*layerd*2
              qstars(:,n_start:n_end,:,:) = 
     +                  RESHAPE(all(is:ie),(/nstars,noccb1,layerd,2/))
              is =ie+1; ie = is - 1 + noccb1 ! bdim
              jsym(n_start:n_end) = NINT(all(is:ie))
              is =ie+1; ie = is - 1 + noccb1 ! bdim
!              IF (ie.ne.id) STOP 'error in mpi_col_dos!'
              ksym(n_start:n_end) = NINT(all(is:ie))

            ELSE ! not l_evp

              is = 6  ; ie = is - 1 + nbdim 
              eig(1:nbdim) = all(is:ie)
              is =ie+1; ie = is - 1 + 4*ntypd*nbdim
              qal(:,:,:) = RESHAPE(all(is:ie),(/4,ntypd,nbdim/))
              is =ie+1; ie = is - 1 + 2*nbdim
              qvac(:,:) = RESHAPE(all(is:ie),(/nbdim,2/))
              is =ie+1; ie = is - 1 + nbdim
              qis(:) = all(is:ie)
              is =ie+1; ie = is - 1 + nbdim*layerd*2
              qvlay(:,:,:) = RESHAPE(all(is:ie),(/nbdim,layerd,2/))
              is =ie+1; ie = is - 1 + nstars*nbdim*layerd*2
           qstars(:,:,:,:)=RESHAPE(all(is:ie),(/nstars,nbdim,layerd,2/))
              is =ie+1; ie = is - 1 + nbdim
              jsym(:) = NINT(all(is:ie))
              is =ie+1; ie = is - 1 + nbdim
              IF (ie.ne.id) STOP 'error in mpi_col_dos!'
              ksym(:) = NINT(all(is:ie))
           
            ENDIF ! not l_evp
c
c--dw   now write k-point data to tmp_dos 
c       in case of evp, write after last communication only
c
            IF ((l_evp.AND.(i == isize-1)).OR.(.not.l_evp)) THEN
              IF ( .not.l_mcd ) THEN
              WRITE (84,rec=ir) bkpt,wk,nbands,eig,
     +                    qal(:,:,:),qvac(:,:),qis(:),
     +                    qvlay(:,:,:),qstars,ksym,jsym
              ELSE
              STOP 'no l_mcd and mpi_col_dos!'
              WRITE (84,rec=ir) bkpt,wk,nbands,eig,
     +                    qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:),
     +                    qstars,ksym,jsym,mcd(:,1:ncored,:)
              ENDIF
            ENDIF
c-new_sl
         ENDDO
      ENDIF
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      DEALLOCATE ( all )
!-----------------------------------------------------------------
      IF (ndir == -3) THEN
      IF (l_evp) STOP 'no eigenvector-parallel & orbcomp'
      id = 5 + neigd*(3 + 2*nsld + 24*natd)
      ALLOCATE ( all(0:id) )
      IF (irank.ne.0) THEN
         all(0)       = nkpt*(jspin-1)+ikpt
         all(1:3)     = bkpt(1:3)
         all(4)       = wk
         all(5)       = real(nbands)
         is = 6  ; ie = is - 1 + neigd
         all(is:ie)   = eig(1:neigd)
         is =ie+1; ie = is - 1 + 2*neigd
         all(is:ie)   = RESHAPE(qvac(:,:),(/neigd*2/))
         is =ie+1; ie = is - 1 + nsld*neigd
         all(is:ie)   = RESHAPE(qintsl(:,:),(/nsld*neigd/))
         is =ie+1; ie = is - 1 + nsld*neigd
         all(is:ie)   = RESHAPE(qmtsl(:,:),(/nsld*neigd/))
         is =ie+1; ie = is - 1 + neigd*23*natd
         all(is:ie)   = RESHAPE(orbcomp(:,:,:),(/neigd*23*natd/))
         is =ie+1; ie = is - 1 + neigd*natd
         IF (ie.ne.id) STOP 'error in mpi_col_dos!'
         all(is:ie)=RESHAPE(qmtp(:,:),(/neigd*natd/))
         CALL MPI_SEND(all,ie,MPI_DOUBLE_PRECISION,0,irank,
     +                                 MPI_COMM_WORLD,ierr)
      ELSE
         DO i = 1,isize-1
           CALL MPI_RECV(all,id,MPI_DOUBLE_PRECISION,i,i,
     +                           MPI_COMM_WORLD,stt,ierr)

            ir        = NINT(all(0))
            bkpt(1:3) = all(1:3)
            wk        = all(4)
            nbands    = NINT(all(5))
            is = 6  ; ie = is - 1 + neigd
            eig(1:neigd) = all(is:ie)
            is =ie+1; ie = is - 1 + 2*neigd
            qvac(:,:) = RESHAPE(all(is:ie),(/neigd,2/))
            is =ie+1; ie = is - 1 + nsld*neigd
            qintsl(:,:) = RESHAPE(all(is:ie),(/nsld,neigd/))
            is =ie+1; ie = is - 1 + nsld*neigd
            qmtsl(:,:) = RESHAPE(all(is:ie),(/nsld,neigd/))
            is =ie+1; ie = is - 1 + neigd*23*natd 
            orbcomp(:,:,:) = RESHAPE(all(is:ie),(/neigd,23,natd/))
            is =ie+1; ie = is - 1 + neigd*natd
            IF (ie.ne.id) STOP 'error in mpi_col_dos!'
            qmtp(:,:)= RESHAPE(all(is:ie),(/neigd,natd/))
c
c--dw   now write k-point data to tmp_dos
c
            IF ( .not.l_mcd ) THEN
              WRITE (129,rec=ir) bkpt,wk,nbands,eig,
     +                          qvac(:,:),qintsl(:,:),qmtsl(:,:),
     +                          orbcomp(:,:,:),qmtp(:,:)
            ELSE
              WRITE (6,FMT=500)
              WRITE (16,FMT=500)
            ENDIF
         ENDDO
      ENDIF

 500  FORMAT('A calculation of an orbital composition of film states
     +       (s,px,py,pz,....) for the l_mcd  case is not implemented.')

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      DEALLOCATE ( all )
      ENDIF ! ndir = -3

      END SUBROUTINE mpi_col_dos
      END MODULE m_mpi_col_dos
