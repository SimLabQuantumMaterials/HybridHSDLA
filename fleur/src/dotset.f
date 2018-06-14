      MODULE m_dotset
c*********************************************************************
c     set up metric matrices for internal coordinates.       m.w.
c*********************************************************************
      CONTAINS
      SUBROUTINE dotset(
     >                  amat,bmat,
     <                  aamat,bbmat)

      IMPLICIT NONE
C     ..
C     .. Array Arguments ..
      REAL, INTENT (IN) ::  amat(3,3), bmat(3,3)
      REAL, INTENT (OUT):: aamat(3,3),bbmat(3,3)
C     ..
C     .. Local Scalars ..
      REAL sa,sb
      INTEGER i,j,k
C     ..
      DO i = 1,3
         DO j = 1,i
            sa = 0.0
            sb = 0.0
            DO k = 1,3
               sa = sa + amat(k,i)*amat(k,j)
               sb = sb + bmat(i,k)*bmat(j,k)
            ENDDO
            aamat(i,j) = sa
            aamat(j,i) = sa
            bbmat(i,j) = sb
            bbmat(j,i) = sb
         ENDDO
      ENDDO

      END SUBROUTINE dotset
      END MODULE m_dotset
