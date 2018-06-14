      MODULE m_mingeselle
      CONTAINS
      SUBROUTINE mingeselle(
     >                      SUB_COMM,n_size,n_rank,nv,
     >                      m_ahelp,matsize,ahelp,
     X                      aa)
c------------------------------------------------------------------+
c                                                                  |
c Transfers the spin-down/spin-up part , upper triangle of the     |
c MT-hamiltonian from the help-array ``ahelp'' to the H-matrix.    |
c For eigenvector-parallelization this needs some communication    |
c between the nodes, since this part is created 'column-wise'      |
c but needed row-wise.                                             |
c                                                                  |
c     n_s(i): number of elements to send to pe #i                  |
c     n_r(i): number of elements to receive from pe #i             |
c     ns_tot,nr_tot : total number of elements to send/receive     |
c     n2_start: pe that has first column of 2nd spin part          |
c     cs_el,cr_el: send and receive elements                       |
c     in_pos(n,i): where to put in the n'th element sent by pe #i  |
c                                                                  |
c------------------------------------------------------------------+
#include"./cpp_double.h"

      IMPLICIT NONE
C ..
C .. Scalar Arguments
      INTEGER, INTENT (IN) :: n_size,n_rank,SUB_COMM
      INTEGER, INTENT (IN) :: m_ahelp,matsize
C ..
C .. Array Arguments
      INTEGER, INTENT (IN) :: nv(2)
      COMPLEX, INTENT (INOUT) :: ahelp(m_ahelp)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: aa(matsize)
#else
      COMPLEX, INTENT (INOUT) :: aa(matsize)
#endif
C ..
C .. Local Scalars
      INTEGER ki,kj,ns_tot,nr_tot,n_p,n2_start,n_help
      INTEGER ns_max,nr_max,n_pos,np_s,np_r,nv_s,ii,i
      INTEGER inext,ifront,req_s,req_r
C ..
C .. Local Arrays
      INTEGER n_s(0:n_size-1),n_send(0:n_size-1)
      INTEGER n_r(0:n_size-1),n_recv(0:n_size-1),ierr(3)
      INTEGER, ALLOCATABLE :: in_pos(:,:)
      COMPLEX, ALLOCATABLE :: cs_el(:,:),cr_el(:),b_b(:),c_help(:,:)

      INCLUDE 'mpif.h'
      INTEGER stt(MPI_STATUS_SIZE)
C ..
c
c kick out the diagonal elements of ahelp
c
      i  = 0
      ii = 0
      DO ki =  n_rank+1, nv(1), n_size
         DO kj = 1,ki - 1
            i  =  i + 1
            ii = ii + 1
            ahelp(i) = ahelp(ii)
         END DO
         ii = ii + 1
      ENDDO
c
c initialize
c
      ns_tot = 0
      nr_tot = 0
      DO n_p = 0,n_size-1
        n_s(n_p) = 0
        n_r(n_p) = 0
      ENDDO
c
c determine number of elements to send to other pe's
c
      n2_start = mod(nv(1),n_size) - 1
      DO ki = 1, nv(1)
        IF ( mod(ki-1,n_size).EQ.n_rank ) THEN
          DO kj = 1, ki-1
            ns_tot = ns_tot + 1
            n_p = mod((kj+n2_start),n_size)
            n_s(n_p) = n_s(n_p) + 1
          ENDDO
        ENDIF
      ENDDO
c
c determine number of elements to receive from other pe's
c
      DO ki = 1, nv(2)
        IF ( mod(ki+nv(1)-1,n_size).EQ.n_rank ) THEN
          DO kj = ki+1, nv(2)
            nr_tot = nr_tot + 1
            n_p = mod(kj-1,n_size)
            n_r(n_p) = n_r(n_p) + 1
          ENDDO
        ENDIF
      ENDDO
c
!      WRITE (*,*) ns_tot,(n_s(n_p),n_p=0,n_size-1)
!      WRITE (*,*) nr_tot,(n_r(n_p),n_p=0,n_size-1)
c
c determine the maximal number of s/r-counts and allocate s/r-arrays
c
      ns_max = 0
      nr_max = 0
      DO n_p = 0,n_size-1
        ns_max = max(ns_max,n_s(n_p))
        nr_max = max(nr_max,n_r(n_p))
      ENDDO
!      WRITE (*,*) ns_max ,nr_max  , n_size, n_rank
      ALLOCATE ( cs_el(ns_max,0:n_size-1),cr_el(nr_max),
     +          in_pos(nr_max,0:n_size-1) )
c
c sort the elements of aahelp-array into the send-arrays
c
      n_help = 0
      DO n_p = 0,n_size-1
        n_send(n_p) = 0
      ENDDO 
      DO ki = 1, nv(1)
        IF ( mod(ki-1,n_size).EQ.n_rank ) THEN
          DO kj = 1, ki-1
            n_help = n_help + 1
            n_p = mod((kj+n2_start),n_size)
            n_send(n_p) = n_send(n_p) + 1
            cs_el(n_send(n_p),n_p) = ahelp(n_help)
          ENDDO
        ENDIF
      ENDDO
      IF (n_help.NE.ns_tot) STOP 'n_help.NE.ns_tot'
      DO n_p = 0,n_size-1
        IF (n_send(n_p).NE.n_s(n_p)) STOP 'n_send.NE.n_s'
      ENDDO
c
c resort send array: rows <-> columns
c
      DO n_p = 0,n_size-1
        nv_s = NINT(sqrt(2.0*n_send(n_p))-0.5)
        ALLOCATE ( c_help(nv_s,nv_s) )

        n_help = 0
        DO ki = 1,nv_s
           DO kj = 1,ki
              n_help = n_help + 1
              c_help(ki,kj) = cs_el(n_help,n_p)
           ENDDO
        ENDDO 
 
        n_help = 0
        DO kj = 1,nv_s
           DO ki = kj ,nv_s
              n_help = n_help + 1
              cs_el(n_help,n_p) = c_help(ki,kj)
           ENDDO
        ENDDO 
 
        DEALLOCATE ( c_help )
      ENDDO
c
c now we look where to put in the received elements
c
      n_pos = 0
      DO n_p = 0,n_size-1
        n_recv(n_p) = 0
      ENDDO
      DO ki = 1, nv(1)+nv(2)
        IF ( mod(ki-1,n_size).EQ.n_rank ) THEN
          DO kj = 1, ki
            n_pos = n_pos + 1 
            IF ( ki.GT.nv(1) ) THEN
              IF ((kj.GT.ki-nv(1)).AND.(kj.LE.nv(1))) THEN
                n_p = mod(kj-1,n_size)
                n_recv(n_p) = n_recv(n_p) + 1
                in_pos(n_recv(n_p),n_p) = n_pos
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO n_p = 0,n_size-1
        IF (n_recv(n_p).NE.n_r(n_p)) STOP 'n_recv.NE.n_s'
      ENDDO
c
c Mandaliet, mandaliet, min geselle kumme niet
c
      ifront = ibefore(n_size,n_rank)
      inext  = iafter (n_size,n_rank)
      DO n_p = 0,n_size-1
c
c determine pe's to send to and to receive from
c
        np_s = mod(inext +n_p,n_size)
        np_r = mod(ifront-n_p,n_size)
        IF (np_r.LT.0) np_r = np_r + n_size
c
c send section: local rows i with mod(i-1,np) = np_s will be sent to proc np_s
c

         IF (np_s.NE.n_rank) THEN
          CALL MPI_ISEND(cs_el(1,np_s),n_send(np_s),
     +                              CPP_MPI_COMPLEX,
     +                              np_s,n_rank,SUB_COMM,req_s,ierr)
!          write (*,*) n_rank,'sends',n_send(np_s),'to',np_s
!          write (*,'(i2,10f10.7)') n_rank,(real(cs_el(ki,np_s)),ki=1,10)
         ENDIF

c
c receive section : local rows i  with mod(i-1,np) = np_r will be received from np_r
c ... skipped, if update matrix from local data:
c
         IF (np_r.NE.n_rank) THEN
           CALL MPI_IRECV(cr_el,n_recv(np_r),CPP_MPI_COMPLEX,
     +                    MPI_ANY_SOURCE,np_r,SUB_COMM,req_r,ierr)
           CALL MPI_WAIT(req_s,stt,ierr)
           CALL MPI_WAIT(req_r,stt,ierr)
!          write (*,*) n_rank,'recvs',ierr,n_p,np_r
!          write(*,*) n_rank,'receives',n_recv(np_r),'from',np_r
!          write (*,'(i2,10f10.7)') n_rank,(real(cr_el(ki)),ki=1,10)
c
c now update the matrix aa()
c
           DO ki=1,n_recv(np_r)
             aa(in_pos(ki,np_r)) = aa(in_pos(ki,np_r)) + cr_el(ki)
           ENDDO
         ELSE
           DO ki=1,n_recv(np_r)
             aa(in_pos(ki,np_r)) = aa(in_pos(ki,np_r)) + cs_el(ki,np_s)
           ENDDO
         ENDIF
!         CALL MPI_BARRIER(SUB_COMM,ierr)
      ENDDO 

      DEALLOCATE (cs_el,cr_el,in_pos)

      END SUBROUTINE mingeselle
!
!-------------------------------------------------------------
!
      INTEGER FUNCTION ibefore(np, p)
!
! Determine (in a ring structure) which is the front process
!
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: np  !  number of processes
      INTEGER, INTENT (IN) :: p   !  current processes

      IF ( p > 0 ) THEN
        ibefore = p-1
      ELSE
        ibefore = np-1
      ENDIF

      END FUNCTION ibefore
!
!-------------------------------------------------------------
!
      INTEGER FUNCTION iafter(np, p)
!
! Determine (in a ring structure) which is the next process
!
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: np  !  number of processes
      INTEGER, INTENT (IN) :: p   !  current processes

      IF ( p < np-1 ) THEN
        iafter = p+1
      ELSE
        iafter = 0
      ENDIF

      END FUNCTION iafter
!
!-------------------------------------------------------------
!
      END MODULE m_mingeselle
