      MODULE m_dosef
c
c--->    obtain dos at ei (here: ef)
c
      CONTAINS
      SUBROUTINE dosef(
     >                 neigd,nkptd,jspd,ntriad,
     >                 ei,nemax,jspins,sfac,ntria,itria,atr,eig)
c
      USE m_trisrt
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd,ntriad
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: ntria
      REAL,    INTENT (IN) :: ei,sfac
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nemax(2)
      INTEGER, INTENT (IN) :: itria(3,ntriad)
      REAL,    INTENT (IN) :: atr(ntriad)
      REAL,    INTENT (IN) :: eig(neigd,nkptd,jspd)
C     ..
C     .. Local Scalars ..
      REAL e1,e2,e21,e3,e31,e32,s
      INTEGER i,jsp,k1,k2,k3,n,neig
      DO 40 jsp = 1,jspins
         neig = nemax(jsp)
         s = 0.0
         DO 30 i = 1,neig
            DO 20 n = 1,ntria
               k1 = itria(1,n)
               k2 = itria(2,n)
               k3 = itria(3,n)
               e1 = eig(i,k1,jsp)
               e2 = eig(i,k2,jsp)
               e3 = eig(i,k3,jsp)
               CALL trisrt(e1,e2,e3,k1,k2,k3)
               IF (e1.LE.-9999.0) GO TO 20
               IF ((ei.LT.e1) .OR. (ei.GE.e3)) GO TO 20
               IF (ei.GT.e2) GO TO 10
c--->    e1<ei<e2
               e31 = e3 - e1
               e21 = e2 - e1
               s = s + 2.*atr(n)* (ei-e1)/ (e31*e21)
               GO TO 20
c--->    e2<ei<e3
   10          e31 = e3 - e1
               e32 = e3 - e2
               s = s + 2.*atr(n)* (e3-ei)/ (e31*e32)
   20       CONTINUE
   30    CONTINUE
!gb         s = (2./jspins)*s
         s = sfac * s
         WRITE (16,FMT=8000) ei,jsp,s
         WRITE (6,FMT=8000) ei,jsp,s
   40 CONTINUE

 8000 FORMAT (/,10x,'density of states at',f12.6,' har for spin',i2,'=',
     +       e20.8,' states/har')

      END SUBROUTINE dosef
      END MODULE m_dosef
