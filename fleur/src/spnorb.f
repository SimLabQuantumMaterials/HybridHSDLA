      MODULE m_spnorb
!*********************************************************************
!     calls soinit to calculate the radial spin-orbit matrix elements:
!     rsopp,rsopdpd,rsoppd,rsopdp
!     and sets up the so - angular matrix elements (soangl)
!     using the functions anglso and sgml.
!*********************************************************************
      CONTAINS
      SUBROUTINE spnorb(
     >                  ntypd,lmaxd,jmtd,jspd,nwdd,nlhd,nlod,
     >                  theta,phi,jspins,ntype,nw,irank,
     >                  jri,lmax,dx,rmsh,el0,ello,nlo,llo,l_dulo,
     >                  ulo_der,vr,spav, 
     <                  rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     <                  rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <                  us,dus,uds,duds,ulos,dulos,uulon,dulon,soangl)

      USE m_soinit 
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,lmaxd,jmtd,jspd,nwdd,nlhd,nlod
      INTEGER, INTENT (IN) :: jspins,ntype,nw,irank
      REAL,    INTENT (IN) :: theta,phi
      LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lmax(ntypd),jri(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntype),llo(nlod,ntype)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
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
      REAL,    INTENT (OUT) ::   us(0:lmaxd,ntype,jspd)
      REAL,    INTENT (OUT) ::  dus(0:lmaxd,ntype,jspd)
      REAL,    INTENT (OUT) ::  uds(0:lmaxd,ntype,jspd)
      REAL,    INTENT (OUT) :: duds(0:lmaxd,ntype,jspd)
      REAL,    INTENT (OUT) ::  ddn(0:lmaxd,ntype,jspd)
      REAL,    INTENT (OUT) ::  ulos(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: dulos(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: uulon(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: dulon(nlod,ntypd,jspd)
      COMPLEX, INTENT (OUT) :: soangl(lmaxd,-lmaxd:lmaxd,2,
     +                                lmaxd,-lmaxd:lmaxd,2)
C     ..
C     .. Local Scalars ..
      INTEGER is1,is2,jspin1,jspin2,l,l1,l2,m1,m2,n
C     ..
C     .. Local Arrays ..
      INTEGER ispjsp(2)
C     ..
C     .. External Functions ..
      REAL sgml
      COMPLEX anglso
      EXTERNAL sgml,anglso
C     ..
      DATA ispjsp/1,-1/

      CALL soinit(
     >            jmtd,jspd,ntypd,nwdd,lmaxd,nlhd,nlod,jri,jspins,
     >            ntype,nw,lmax,dx,rmsh,el0,ello,nlo,llo,l_dulo,
     >            ulo_der,vr,spav,
     <            rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     <            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <            us,dus,uds,duds,ulos,dulos,uulon,dulon)
c
      IF (irank.EQ.0) THEN
        DO n = 1,ntype
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,1),l=1,3)
         WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,1),l=1,3)
         WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,1),l=1,3)
         WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,1),l=1,3)
        ENDDO
      ENDIF
 8000 FORMAT (' spin - orbit parameter HR  ')
 8001 FORMAT (8f8.4)
 9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')
c

      IF ((abs(theta).LT.0.00001).AND.(abs(phi).LT.0.00001)) THEN
c
c       TEST for real function sgml(l1,m1,is1,l2,m2,is2)
c
        DO l1 = 1,lmaxd
          DO l2 = 1,lmaxd
            DO jspin1 = 1,2
              DO jspin2 = 1,2
                is1=ispjsp(jspin1)
                is2=ispjsp(jspin2)
                DO m1 = -l1,l1,1
                  DO m2 = -l2,l2,1
                    soangl(l1,m1,jspin1,l2,m2,jspin2) =
     +                    cmplx(sgml(l1,m1,is1,l2,m2,is2),0.0)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
   
      ELSE
c
c       TEST for complex function anglso(teta,phi,l1,m1,is1,l2,m2,is2)
c 
        DO l1 = 1,lmaxd
          DO l2 = 1,lmaxd
            DO jspin1 = 1,2
              DO jspin2 = 1,2
                is1=ispjsp(jspin1)
                is2=ispjsp(jspin2)
c
                DO m1 = -l1,l1,1
                  DO m2 = -l2,l2,1
                    soangl(l1,m1,jspin1,l2,m2,jspin2) =
     =                    anglso(theta,phi,l1,m1,is1,l2,m2,is2)
                  ENDDO
                ENDDO
c
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c
      ENDIF

      IF (irank.EQ.0) THEN
        WRITE (6,FMT=8002)
        DO jspin1 = 1,2
          DO jspin2 = 1,2
            WRITE (6,FMT=*) 'd-states:is1=',jspin1,',is2=',jspin2
            WRITE (6,FMT='(7x,7i8)') (m1,m1=-3,3,1)
            WRITE (6,FMT=8003) (m2, (soangl(3,m1,jspin1,3,m2,jspin2),
     +                                           m1=-3,3,1),m2=-3,3,1)
          ENDDO
        ENDDO
      ENDIF

 8002 FORMAT (' so - angular matrix elements ')
 8003 FORMAT (i8,14f8.4)

      END SUBROUTINE spnorb
      END MODULE m_spnorb
