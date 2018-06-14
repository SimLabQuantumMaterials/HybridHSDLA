      MODULE m_mpi_col_eig
!
! collects the eigenvalue-information necessary for the determination 
! of the fermi energy from the nodes and writes it to the eig-file of
! pe 0; only for parallel applications without direct global access
! files.                                                        gb 01
!     
      CONTAINS
      SUBROUTINE mpi_col_eig(
     >                       irank,isize,jspd,neigd,lmaxd,ntypd,nlod,
     >                  l_ss,l_noco,l_J,jspin,form66,eonly,film,nkpt,
     >                       nrec,el,evac,ello,bkpt,wtkpt,ne,nv,
     >                       nmat,eig)

      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT (IN) :: jspd,neigd,lmaxd,ntypd,nlod
      INTEGER, INTENT (IN) :: irank,isize,jspin,nkpt
      LOGICAL, INTENT (IN) :: l_ss,l_noco,form66,eonly,film,l_J
      INTEGER, INTENT (IN) :: nrec,ne,nmat
      INTEGER, INTENT (IN) :: nv(jspd)
      REAL,    INTENT (IN) :: eig(neigd),bkpt(3)
      REAL,    INTENT (IN) :: el(0:lmaxd,ntypd,jspd),wtkpt
      REAL,    INTENT (IN) :: ello(nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: evac(2,jspd)


      INTEGER i,nb,nrec_l,ne_l,l,n,req_s,req_r
      INTEGER ierr(3),stt(MPI_STATUS_SIZE)
      REAL evac_sv
      REAL, ALLOCATABLE :: all(:)

      ALLOCATE ( all(neigd+7) )
      IF (irank.ne.0) THEN
         all(1)       = real(nrec)
         all(2)       = real(ne)
         all(3)       = wtkpt
         all(4:6)     = bkpt(1:3)
         all(7:neigd+6) = eig(1:neigd)
         all(neigd+7) = evac(1,1)
         CALL MPI_ISEND(all,neigd+7,MPI_DOUBLE_PRECISION,0,irank,
     +                             MPI_COMM_WORLD,req_s,ierr)
         CALL MPI_WAIT(req_s,stt,ierr)
      ELSE
         evac_sv = evac(1,1) 
         DO i = 1,isize-1
           CALL MPI_IRECV(all,neigd+7,MPI_DOUBLE_PRECISION,i,i,
     +                               MPI_COMM_WORLD,req_r,ierr)
           CALL MPI_WAIT(req_r,stt,ierr)
           evac(1,1) = all(neigd+7)
           nrec_l = int( all(1) )
           IF (.not.l_J) write(*,*) 'now',nrec+i,' instead of',nrec_l
           nrec_l = nrec + i
           ne_l =   int( all(2) )
           IF (.not.l_J) WRITE (6,FMT=8010) ne_l,i
           IF (.not.l_J) WRITE (6,FMT=8020) (all(6+nb),nb=1,ne_l)
!
! --> write info on eig-file of pe 0
!
           IF (form66) THEN
             IF (.NOT.eonly) THEN
               WRITE (66,FMT=8110) ((el(l,n,jspin),l=0,lmaxd),n=1,ntypd)
               IF (film) WRITE (66,FMT=8110) (evac(l,jspin),l=1,2)
             ENDIF
             WRITE (66,FMT=8070) all(4:6),all(3),ne_l,nv,nmat
             WRITE (66,FMT=8080) (all(nb),nb=7,6+ne_l)
             IF (.NOT.eonly) THEN
               WRITE (66,*) 'Sorry, no formt. EVs on parallel machines'
!              WRITE (66,FMT=8090) (k1(k,jspin),k2(k,jspin),k3(k,jspin)
!    +                             ,k=1,nv(jspin))
!              DO j = 1,ne
!                WRITE (66,FMT=8100) (z(k,j),k=1,nv(jspin))
!              ENDDO
             ENDIF
 8110        FORMAT (1p,4e20.13)
 8070        FORMAT (3f10.6,f12.8,4i5)
 8080        FORMAT (6f15.6)
 8090        FORMAT (12i6)
 8100        FORMAT (1p,4e20.13)

           ELSE
             IF (l_J) THEN
               WRITE (66,rec=nrec_l) 
     +                      all(4:6),all(3),ne_l,all(7:neigd+6)
             ELSE
               IF (l_ss) THEN
                 WRITE (66,rec=nrec_l)
     +              el,evac,ello,all(4:6),all(3),ne_l,nv,nmat,
     +              all(7:neigd+6)
               ELSEIF (l_noco) THEN
                 WRITE (66,rec=nrec_l)
     +              el,evac,ello,all(4:6),all(3),ne_l,nv(jspin),nmat,
     +              all(7:neigd+6)
               ELSE
                 WRITE (66,rec=nrec_l) el(:,:,jspin),evac(:,jspin),
     +              ello(:,:,jspin),all(4:6),all(3),ne_l,
     +              nv(jspin),nmat,all(7:neigd+6)
              ENDIF
             ENDIF
           ENDIF

           IF (mod(nrec_l,nkpt)==0) GOTO 100
         ENDDO
         evac(1,1) = evac_sv
      ENDIF

      !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  100 DEALLOCATE ( all )
 8010 FORMAT (' the',i4,' eigenvalues on node',i4,' are:')
 8020 FORMAT (5x,5f12.6)

      END SUBROUTINE mpi_col_eig
      END MODULE m_mpi_col_eig
