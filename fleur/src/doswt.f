      MODULE m_doswt
c
c     calculates the weights for each k-point for integrating functions
c     of k.  the array w has beeen cleared before entering.
c
      CONTAINS
      SUBROUTINE doswt(
     >                 neigd,nkptd,jspd,ntriad,
     >                 ei,nemax,jspins,ntria,itria,atr,eig,
     <                 w)
      USE m_trisrt
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd,ntriad
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: ntria
      REAL,    INTENT (IN) :: ei
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nemax(2)
      INTEGER, INTENT (IN) :: itria(3,ntriad)
      REAL,    INTENT (IN) :: atr(ntriad)
      REAL,    INTENT (IN) :: eig(neigd,nkptd,jspd)
      REAL,    INTENT (OUT):: w(neigd,nkptd,jspd)
C     ..
C     .. Local Scalars ..
      INTEGER jsp,i,n
      INTEGER k1,k2,k3
      INTEGER neig
      REAL    e1,e2,e3
      REAl    ee,e32,e31,e21,s

      DO 140 jsp = 1,jspins
         neig = nemax(jsp)
         DO 130 i = 1,neig
            DO 120 n = 1,ntria
               k1 = itria(1,n)
               k2 = itria(2,n)
               k3 = itria(3,n)
               e1 = eig(i,k1,jsp)
               e2 = eig(i,k2,jsp)
               e3 = eig(i,k3,jsp)
               CALL trisrt(e1,e2,e3,k1,k2,k3)
               IF (e1.LE.-9999.0) GO TO 120
               IF (ei.LE.e1) GO TO 120
               IF (ei.GE.e3) GO TO 110
               IF (ei.GT.e2) GO TO 100
c--->    e1<ei<e2
               ee = ei - e1
               e31 = ee/ (e3-e1)
               e21 = ee/ (e2-e1)
               s = atr(n)*e31*e21/3.0
               w(i,k1,jsp) = w(i,k1,jsp) + s* (3.0-e21-e31)
               w(i,k2,jsp) = w(i,k2,jsp) + s*e21
               w(i,k3,jsp) = w(i,k3,jsp) + s*e31
               GO TO 120
c--->    e2<ei<e3
  100          ee = e3 - ei
               e31 = ee/ (e3-e1)
               e32 = ee/ (e3-e2)
               s = atr(n)/3.
               w(i,k1,jsp) = w(i,k1,jsp) + s* (1.-e31*e31*e32)
               w(i,k2,jsp) = w(i,k2,jsp) + s* (1.-e31*e32*e32)
               w(i,k3,jsp) = w(i,k3,jsp) + s* (1.-e31*e32* (3.-e31-e32))
               GO TO 120
c--->    e3<ei
  110          s = atr(n)/3.
               w(i,k1,jsp) = w(i,k1,jsp) + s
               w(i,k2,jsp) = w(i,k2,jsp) + s
               w(i,k3,jsp) = w(i,k3,jsp) + s
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE

      END SUBROUTINE doswt
      END MODULE m_doswt
