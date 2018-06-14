      MODULE m_dotir

      IMPLICIT NONE
      REAL,    PRIVATE :: s,x
      INTEGER, PRIVATE :: i,j

      CONTAINS
c*********************************************************************
      REAL FUNCTION dotirl(f,g,aamat)
c     calculates the dot product of real space vectors f,g given in
c     internal coordinates.                            m.w.
c*********************************************************************
C     ..
C     .. Array Arguments ..
      REAL, INTENT (IN) :: f(3),g(3),aamat(3,3)
C     ..
      s = 0.0
      DO i = 1,3
         x = 0.0
         DO j = 1,3
            x = x + aamat(j,i)*g(j)
         ENDDO
         s = s + f(i)*x
      ENDDO
      dotirl = s
      END FUNCTION dotirl
c*********************************************************************
      REAL FUNCTION dotirp(f,g,bbmat)
c     calculated the dot product of reciprocalspace vectors f,g given
c     in internal coordinates.                        m.w.
c*********************************************************************
C     ..
C     .. Array Arguments ..
      REAL, INTENT(IN) :: f(3),g(3),bbmat(3,3)
C     ..
      s = 0.0
      DO i = 1,3
         x = 0.0
         DO j = 1,3
            x = x + bbmat(j,i)*g(j)
         ENDDO
         s = s + f(i)*x
      ENDDO
      dotirp = s
      END FUNCTION dotirp

      END MODULE m_dotir
