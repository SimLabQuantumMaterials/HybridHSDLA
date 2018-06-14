      MODULE m_mpi_col_eigJ
!
! collects the eigenvalue-information necessary for the determination 
! of the fermi energy from the nodes and writes it to the eig-file of
! pe 0; only for parallel applications without direct global access
! files.                                                        gb 01
!     
      CONTAINS
      SUBROUTINE mpi_col_eigJ(
     >                       irank,isize,nkptd,neigd,nkpt,nkpt1,eig_l,
     <                       bk,wtkpt,ne,eig)

      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT (IN) :: irank,isize,nkptd,neigd,nkpt,nkpt1

      REAL,    INTENT (IN)  :: eig_l(neigd+5,nkpt1)
      REAL,    INTENT (OUT) :: eig(neigd,nkptd)
      REAL,    INTENT (OUT) :: wtkpt(nkptd),bk(3,nkptd)
      INTEGER, INTENT (OUT) :: ne(nkptd)


      INTEGER i,nb,ne_l,l,n,req_s,req_r,nkpt_l,ik,ik_l
      INTEGER ierr(3),stt(MPI_STATUS_SIZE)
      REAL evac_sv
      REAL, ALLOCATABLE :: all(:,:)

      nkpt_l = CEILING( (nkpt - (irank + 0.999999) )/ isize)
      IF (irank.ne.0) THEN
         n = nkpt_l * (neigd+5)
         CALL MPI_ISEND(eig_l,n,MPI_DOUBLE_PRECISION,0,irank,
     +                             MPI_COMM_WORLD,req_s,ierr)
         CALL MPI_WAIT(req_s,stt,ierr)
      ELSE
         DO ik = 1, nkpt_l
           bk(1:3,ik) = eig_l(1:3,ik)
           wtkpt(ik)  = eig_l(4,ik)
           ne(ik)     = NINT(eig_l(5,ik))
           eig(1:neigd,ik) = eig_l(6:neigd+5,ik)
         ENDDO
         ik = nkpt_l
         DO i = 1,isize-1
           nkpt_l = CEILING( (nkpt - (i + 0.999999) )/ isize)
           n = nkpt_l * (neigd+5)
           ALLOCATE ( all(neigd+5,nkpt_l) )
           CALL MPI_IRECV(all,n,MPI_DOUBLE_PRECISION,i,i,
     +                               MPI_COMM_WORLD,req_r,ierr)
           CALL MPI_WAIT(req_r,stt,ierr)

           DO ik_l = 1, nkpt_l
             ik = ik + 1
             bk(1:3,ik) = all(1:3,ik_l)
             wtkpt(ik)  = all(4,ik_l)
             ne(ik)     = NINT(all(5,ik_l))
             eig(1:neigd,ik) = all(6:neigd+5,ik_l)
           ENDDO
           DEALLOCATE (all)
         ENDDO
      ENDIF

      END SUBROUTINE mpi_col_eigJ
      END MODULE m_mpi_col_eigJ
