      MODULE m_dosint
c
c     integrated dos to ei
c
      CONTAINS
      SUBROUTINE dosint(
     >                  neigd,nkptd,jspd,ntriad,
     >                  ei,nemax,jspins,sfac,ntria,itria,atr,eig,
     <                  ct)
      USE m_trisrt
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd,ntriad
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: ntria
      REAL,    INTENT (IN) :: ei,sfac
      REAL,    INTENT (OUT):: ct
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nemax(2)
      INTEGER, INTENT (IN) :: itria(3,ntriad)
      REAL,    INTENT (IN) :: atr(ntriad)
      REAL,    INTENT (IN) :: eig(neigd,nkptd,jspd)
C     ..
C     .. Local Scalars ..
      INTEGER jsp,i,n
      INTEGER k1,k2,k3
      INTEGER neig
      REAL    e1,e2,e3
      REAl    ee,e32,e31,e21,s
c
      s = 0.0
      DO 90 jsp = 1,jspins
         neig = nemax(jsp)
         DO 80 i = 1,neig
            DO 70 n = 1,ntria
               k1 = itria(1,n)
               k2 = itria(2,n)
               k3 = itria(3,n)
               e1 = eig(i,k1,jsp)
               e2 = eig(i,k2,jsp)
               e3 = eig(i,k3,jsp)
c     write(16,*) 'eig=',i,', ntria=',n,', e1,e2,e3=',e1,e2,e3
               CALL trisrt(e1,e2,e3,k1,k2,k3)
               IF (e1.LE.-9999.0) GO TO 70
               IF (ei.LE.e1) GO TO 70
               IF (ei.GE.e3) GO TO 60
               IF (ei.GT.e2) GO TO 50
c--->    e1<ei<e2
               e21 = e2 - e1
               e31 = e3 - e1
               ee = ei - e1
               s = s + atr(n)*ee*ee/ (e21*e31)
               GO TO 70
c--->    e2<ei<e3
   50          e31 = e3 - e1
               e32 = e3 - e2
               ee = e3 - ei
               s = s + atr(n)* (1.-ee*ee/ (e31*e32))
               GO TO 70
c--->    e3<ei
   60          s = s + atr(n)
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
cjr      ct=2.*s
!gb      ct = (2./jspins)*s
      ct = sfac * s

      END SUBROUTINE dosint
      END MODULE m_dosint
