      MODULE m_ccdnup
c     *******************************************************
c     *****   set up the core densities for compounds.  *****
c     *****   in accordanse to d.d.koelling's cored     *****
c     *******************************************************
      CONTAINS
      SUBROUTINE ccdnup(
     >                  jmtd,jspd,msh,nlhd,ntypd,jspins,jatom,
     >                  jri,neq,dx,rho,rmsh,rmt,
     >                  sume,vrs,rhochr,rhospn,
     <                  tecs,qints)

      USE m_intgr, ONLY : intgr3
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,jspd,msh,nlhd,ntypd
      INTEGER, INTENT (IN) :: jatom,jspins
      REAL,    INTENT (IN) :: sume
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),neq(ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd),rmt(ntypd)
      REAL,    INTENT (IN) :: rhochr(msh),rhospn(msh)
      REAL,    INTENT (IN) :: vrs(jmtd,ntypd,jspd)
      REAL,    INTENT (OUT) :: tecs(ntypd,jspd)
      REAL,    INTENT (OUT) :: qints(ntypd,jspd)
      REAL,    INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
C     ..
C     .. Local Scalars ..
      REAL d,dxx,q,rad,rhs,sfp
      INTEGER i,j,jspin,nm,nm1
C     ..
C     .. Local Arrays ..
      REAL rhoc(msh),rhoss(msh)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp
C     ..
      sfp = 4. * sqrt( atan (1.) )
      nm = jri(jatom)
c     ---->update spherical charge density rho with the core density.
c     ---->for spin-polarized (jspins=2), take only half the density
      DO jspin = 1,jspins
         IF (jspins.EQ.2 .AND. jspin.EQ.1) THEN
            DO j = 1,msh
               rhoss(j) = rhochr(j) - rhospn(j)
            END DO
         ELSE IF (jspins.EQ.2 .AND. jspin.EQ.2) THEN
            DO j = 1,msh
               rhoss(j) = rhochr(j) + rhospn(j)
            END DO
c jspins=1
         ELSE
            DO j = 1,msh
               rhoss(j) = rhochr(j)
            END DO
c
         END IF
c
         DO 10 j = 1,nm
            rhoc(j) = rhoss(j)/jspins
            rho(j,0,jatom,jspin) = rho(j,0,jatom,jspin) + rhoc(j)/sfp
   10    CONTINUE
         DO 30 i = 1,nm
            rhoc(i) = rhoc(i)*vrs(i,jatom,jspin)/rmsh(i,jatom)
   30    CONTINUE
         CALL intgr3(rhoc,rmsh(1,jatom),dx(jatom),nm,rhs)
         tecs(jatom,jspin) = sume/jspins - rhs
         WRITE (6,FMT=8010) jatom,jspin,tecs(jatom,jspin),sume/jspins
         WRITE (16,FMT=8010) jatom,jspin,tecs(jatom,jspin),sume/jspins
c     write(17) tec
   40    CONTINUE
c     ---> simpson integration
         dxx = dx(jatom)
         d = exp(dx(jatom))
         rad = rmt(jatom)
         q = rad*rhoss(nm)/2.
         DO 50 nm1 = nm + 1,msh - 1,2
            rad = d*rad
            q = q + 2*rad*rhoss(nm1)
            rad = d*rad
            q = q + rad*rhoss(nm1+1)
   50    CONTINUE
         q = 2*q*dxx/3
c+sb
         WRITE (6,FMT=8000) q/jspins
         WRITE (16,FMT=8000) q/jspins
c-sb
         qints(jatom,jspin) = q*neq(jatom)

      END DO ! end-do-loop jspins

 8000 FORMAT (f20.8,' electrons lost from core.')
 8010 FORMAT (10x,'atom type',i3,'  (spin',i2,') ',/,10x,
     +       'kinetic energy=',e20.12,5x,'sum of the eigenvalues=',
     +       e20.12)

      END SUBROUTINE ccdnup
      END MODULE m_ccdnup
