      SUBROUTINE fitchk(f1,f2,n,av,rms,dmx)
c     ************************************************
c     compare functions f1 and f2
c     ************************************************
C     .. Scalar Arguments ..

      IMPLICIT NONE
      REAL av,dmx,rms
      INTEGER n
C     ..
C     .. Array Arguments ..
      REAL f1(n),f2(n)
C     ..
C     .. Local Scalars ..
      REAL d
      INTEGER i
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,sqrt
C     ..
      av = 0.
      rms = 0.
      dmx = 0.
      DO 10 i = 1,n
         av = av + f1(i)
         d = (f1(i)-f2(i))**2
         dmx = max(d,dmx)
         rms = rms + d
   10 CONTINUE
      av = av/n
      IF (abs(av).LT.1.e-30) THEN
         rms = 0.
         dmx = 0.
         RETURN
      END IF
      rms = sqrt(rms/n)/av*100.
      dmx = sqrt(dmx)/av*100.
      RETURN
      END
