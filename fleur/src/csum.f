      COMPLEX FUNCTION csum(n,cvec,incx)

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER incx,n
C     ..
C     .. Array Arguments ..
      COMPLEX cvec(n)
C     ..
C     .. Local Scalars ..
      INTEGER i
C     ..
      csum = cmplx(0.0,0.0)
      DO 10 i = 1,n,incx
         csum = csum + cvec(i)
   10 CONTINUE
      RETURN
      END
