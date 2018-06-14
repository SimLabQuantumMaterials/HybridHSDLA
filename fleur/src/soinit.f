      MODULE m_soinit
c
c**********************************************************************
c     1. generates radial spin-orbit matrix elements:sorad
c     2. generates spin-angular spin-orbit matrix   :soorb (not implemented)
c**********************************************************************
c
      CONTAINS
      SUBROUTINE soinit(
     >                  jmtd,jspd,ntypd,nwdd,lmaxd,nlhd,nlod,jri,jspins,
     >                  ntype,nw,lmax,dx,rmsh,el0,ello,nlo,llo,l_dulo,
     >                  ulo_der,vr,spav,
     <                  rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     <                  rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <                  us,dus,uds,duds,ulos,dulos,uulon,dulon)

      USE m_sorad
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,jspd,ntypd,nwdd,lmaxd,nlod
      INTEGER, INTENT (IN) :: jspins,ntype,nw,nlhd
      LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lmax(ntypd),jri(ntypd),ulo_der(nlod,ntypd)
      INTEGER, INTENT (IN) :: nlo(ntype),llo(nlod,ntype)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (IN) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,jspd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      REAL,    INTENT (OUT) :: rsopp  (ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsoppd (ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsopdp (ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsopdpd(ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsoplop (ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsoplopd(ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsopdplo(ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsopplo (ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsoploplop(ntypd,nlod,nlod,2,2)
      REAL,    INTENT (OUT) ::   us(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  dus(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  uds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) :: duds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  ddn(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  ulos(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: dulos(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: uulon(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: dulon(nlod,ntypd,jspd)
C     ..
C     .. Local Scalars ..
      INTEGER i,jspin,n
C     ..
C     .. Local Arrays ..
      REAL vr0(jmtd,jspd)
C     ..
      DO n = 1,ntype
c
          DO jspin = 1,jspins
c
              DO i = 1,jri(n)
                  vr0(i,jspin) = vr(i,0,n,jspin)
              END DO
              DO i = jri(n) + 1,jmtd
                  vr0(i,jspin) = 0.0
              END DO
c
          END DO
c
          CALL sorad(
     >      jmtd,jspd,ntypd,nwdd,lmaxd,nlod,nlo(n),llo(1,n),
     >      l_dulo(1,n),ulo_der(1,n),ello,
     >      jspins,n,nw,vr0,rmsh(1,n),dx(n),jri(n),lmax(n),el0,spav,
     <      rsopp,rsopdpd,rsoppd,rsopdp,ddn,
     <      rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <      us,dus,uds,duds,ulos,dulos,uulon,dulon)
c
      END DO ! end-do-loop : ntype
 
      END SUBROUTINE soinit
      END MODULE m_soinit
