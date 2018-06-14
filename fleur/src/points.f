      SUBROUTINE points(x,n)
c     *********************************************************
c     generate random points, in internal coordinates,
c     within the unit cell omega-tilda
c     *********************************************************

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER n
C     ..
C     .. Array Arguments ..
      REAL x(3,n)
C     ..
C     .. Local Scalars ..
      REAL r
      INTEGER i,j
C     ..
C     .. External Functions ..
      REAL qranf
      EXTERNAL qranf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
      r = sqrt(13.)
      j = 1
      DO 10 i = 1,n
         x(1,i) = qranf(r,j)
         x(2,i) = qranf(r,j)
         x(3,i) = qranf(r,j) - 0.5
   10 CONTINUE
      RETURN
      END
