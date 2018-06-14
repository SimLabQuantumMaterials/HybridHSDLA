      MODULE m_sointg
c*********************************************************************
c     compute radial spin-orbit integrant
c*********************************************************************
      CONTAINS
      SUBROUTINE sointg(
     >                  e,vr,v0,r0,dx,jri,jmtd,c,jspins,
     <                  vso)
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jri,jspins,jmtd
      REAL,    INTENT (IN) :: dx,e,r0,c
C     ..
C     .. Array Arguments ..
      REAL, INTENT (IN)  ::  v0(jri)
      REAL, INTENT (IN)  ::  vr(jmtd,jspins)
      REAL, INTENT (OUT) :: vso(jmtd,2)
C     ..
C     .. Local Scalars ..
      REAL dr,dxv,r,cin2
      INTEGER i,jspin
C     ..
C     .. Local Arrays ..
      REAL dv(jri),xmrel(jri)
C     ..
C     .. External Subroutines ..
      EXTERNAL diff3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp
C     ..
      cin2 = 1.0/c**2
      dr = exp(dx)
      r = r0
c
c---> relativistic mass:
      DO i = 1,jri
          xmrel(i) = 1. + (e-v0(i)/r)/2.*cin2
          r = r*dr
      END DO
c---> potential derivative (on logarithmic mesh) : v0 := r*v


      CALL diff3(
     >           vr(1,1),jri,dx,
     <           dv)
c
      r = r0
      DO i = 1,jri
          dxv = (dv(i) - vr(i,1))/(r**3)
          vso(i,2) = cin2*dxv/(4.*xmrel(i)**2)
          r = r*dr
      END DO

      IF (jspins.EQ.2) THEN
        CALL diff3(
     >             vr(1,jspins),jri,dx,
     <             dv)
c
        r = r0
        DO i = 1,jri
            dxv = (dv(i) - vr(i,jspins))/(r**3)
            vso(i,1) = cin2*dxv/(4.*xmrel(i)**2)
            r = r*dr
        END DO
      ELSE 
        DO i = 1,jri
            vso(i,1) =  vso(i,2)
        END DO
      ENDIF

      END SUBROUTINE sointg
      END MODULE m_sointg
