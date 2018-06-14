      MODULE m_stmix
c
c      straight mixing,     r.pentcheva, iff, 1996
c
c     sm   : input charge density of iteration m
c     sm1  : input charge density of iteration m+1
c     fsm  : output minus input charge densityof iteration m
c
      CONTAINS
      SUBROUTINE stmix(
     >                 mmap,jspd,n_u,
     >                 nmap,nmaph,jspins,l_noco,alpha,spinf,fsm,
     =                 sm)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      LOGICAL, INTENT (IN) :: l_noco
      INTEGER, INTENT (IN) :: jspins,nmaph,nmap,mmap,jspd,n_u
      REAL,    INTENT (IN) :: alpha,spinf
C     ..
C     .. Array Arguments ..
      REAL fsm(mmap),sm(mmap)
C     ..
C     .. Local Scalars ..
      INTEGER imap
      REAL tol_6
C     ..
C     .. Data statements ..
      DATA tol_6 /1.0e-6/
C     ..
c
      WRITE (16,FMT='(a)') 'STRAIGHT MIXING'
      IF (jspins.EQ.1) WRITE (16,FMT='(a,2f10.5)')
     +    'charge density mixing parameter:',alpha
      IF (jspins.EQ.2) WRITE (16,FMT='(a,2f10.5)')
     +    'spin density mixing parameter:',alpha*spinf
      IF ( abs(spinf-1.0e0).LE.tol_6 .OR. jspins.eq.1 ) then
c     --> perform simple mixing 
c
c        sm1 = sm + alpha * F(sm)

         DO imap = 1,nmap
            sm(imap) = sm(imap) + alpha*fsm(imap)
         END DO
         RETURN
      ELSE
c     -->perform simple mixing with the mixing parameters
c        for charge and spin
c
c       sm1+/_ = (sm+/_) + alpha* F(sm)
c                +/-0.5alpha(spinf-1)( F(sm+) + F(sm-) )

         DO imap = 1,nmaph
            sm(imap) = sm(imap) + alpha*fsm(imap) 
     +            + alpha/2.0*(spinf-1.0)*(fsm(imap) - fsm(imap+nmaph))
         ENDDO

         DO imap = nmaph+1,2*nmaph
            sm(imap) = sm(imap) + alpha*fsm(imap) 
     +            + alpha/2.0*(spinf-1.0)*(fsm(imap) - fsm(imap-nmaph))
         ENDDO
         IF (l_noco) THEN
            DO imap = 2*nmaph+1, nmap - 98*jspins*n_u
               sm(imap) = sm(imap) + alpha*spinf*fsm(imap) 
            ENDDO
         ENDIF
         IF ( n_u > 0 )  THEN
           DO imap = nmap - 98*jspins*n_u + 1, nmap - 98*n_u 
            sm(imap) = sm(imap) + alpha*fsm(imap) 
     +            + alpha/2.0*(spinf-1.0)*(fsm(imap) - fsm(imap+98*n_u))
           ENDDO
           DO imap = nmap - 98*n_u + 1, nmap
            sm(imap) = sm(imap) + alpha*fsm(imap) 
     +            + alpha/2.0*(spinf-1.0)*(fsm(imap) - fsm(imap-98*n_u))
           ENDDO
         ENDIF
      END IF
          
      END SUBROUTINE stmix
      END MODULE m_stmix
