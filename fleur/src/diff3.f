      SUBROUTINE diff3(
     >                 f,jri,dx,
     <                 df)
c
c********************************************************************
c     call       diff3(v0,jri,dx,dv)
c...............................................................diff3
c     differetiation via 3-points
c********************************************************************
c
      IMPLICIT NONE
C     .. 
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jri
      REAL,    INTENT (IN) :: dx
C     ..
C     .. Array Arguments ..
      REAL, INTENT (IN)  ::  f(jri)
      REAL, INTENT (OUT) :: df(jri)
C     ..
C     .. Local Scalars ..
      INTEGER i
      REAL tdx_i
C     ..
      tdx_i = 1./(2.*dx)
c
c---> first point
      df(1) = -tdx_i * (-3.*f(1)+4.*f(2)-f(3))
c
c---> central point formula in charge
      DO i = 2,jri - 1
          df(i) = tdx_i * (f(i+1)-f(i-1))
      END DO
c
c---> last point
      df(jri) = tdx_i * (3.*f(jri)-4.*f(jri-1)+f(jri-2))
c
      RETURN
      END
