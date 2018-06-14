      REAL FUNCTION ssum(n,x,incx)

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER incx,n
C     ..
C     .. Array Arguments ..
      REAL x(n)
C     ..
C     .. Local Scalars ..
      INTEGER i
C     ..
      ssum = 0.
      DO 10 i = 1,n,incx
         ssum = ssum + x(i)
   10 CONTINUE
      RETURN
      END
