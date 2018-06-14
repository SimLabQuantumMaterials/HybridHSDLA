      SUBROUTINE rad_ovlp(
     >           jspd,ntypd,lmaxd,nlhd,jmtd,
     >           jspins,ntype,lmax,jri,rmsh,dx,vr,epar,
     <           uun21,udn21,dun21,ddn21)
c***********************************************************************
c calculates the overlap of the radial basis functions with different
c spin directions. These overlapp integrals are needed to calculate
c the contribution to the hamiltonian from the constant constraint
c B-field.
c
c Philipp Kurz 2000-04
c***********************************************************************

      USE m_int21
      USE m_radfun
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,ntypd,lmaxd,nlhd,jmtd
      INTEGER, INTENT (IN) :: jspins,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT  (IN):: lmax(ntypd),jri(ntypd)
      REAL,    INTENT  (IN):: rmsh(jmtd,ntypd),dx(ntypd)
      REAL,    INTENT  (IN):: epar(0:lmaxd,ntypd,jspd)
      REAL,    INTENT  (IN):: vr(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (OUT):: uun21(0:lmaxd,ntypd),udn21(0:lmaxd,ntypd)
      REAL,    INTENT (OUT):: dun21(0:lmaxd,ntypd),ddn21(0:lmaxd,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER itype,l,ispin,noded,nodeu
      REAL ddn,duds,dus,uds,us,wronk
C     ..
C     .. Local Arrays ..
      REAL f(jmtd,2,0:lmaxd,jspd),g(jmtd,2,0:lmaxd,jspd)
C     ..
      DO itype = 1,ntype
         DO l = 0,lmax(itype)
            DO ispin = 1,jspins
               CALL radfun(
     >              l,epar(l,itype,ispin),vr(1,0,itype,ispin),
     >              jri(itype),rmsh(1,itype),dx(itype),jmtd,
     <              f(1,1,l,ispin),g(1,1,l,ispin),us,dus,uds,duds,ddn,
     <              nodeu,noded,wronk)
            ENDDO
            CALL int_21(
     >                  f,g,rmsh(1,itype),dx(itype),jri(itype),
     >                  jmtd,lmaxd,jspd,l,
     <                  uun21(l,itype),udn21(l,itype),
     <                  dun21(l,itype),ddn21(l,itype))
         ENDDO
      ENDDO

      END

