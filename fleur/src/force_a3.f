      MODULE m_force_a3
      CONTAINS
      SUBROUTINE force_a3(jmtd,memd,nlhd,ntypd,jspd,ntypsd,natd,
     >                    jri,mlh,nlh,nmem,ntype,jspins,ntypsy,
     >                    zatom,dx,rmt,rmsh,clnu,llh,rho,vr,neq,
     >                    l_geo,
     X                    force)
c ************************************************************
c Hellman-Feynman force contribution a la Rici et al.
c ************************************************************
c
      USE m_intgr, ONLY : intgr3
      USE m_constants, ONLY : pimach
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,memd,nlhd,ntypd,jspd,ntypsd,natd
      INTEGER, INTENT (IN) :: ntype,jspins
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nlh(ntypsd),mlh(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: zatom(ntypd),dx(ntypd),rmt(ntypd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) ::  vr(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
C     ..
C     .. Local Scalars ..
      COMPLEX czero
      REAL a3_1,a3_2,s34,s38,w,pi
      INTEGER i,ir,jsp,lh,lindex,mem,mindex,n,nd,na
C     ..
C     .. Local Arrays ..
      COMPLEX forc_a3(3),grd1(3,-1:1),gv(3)
      REAL rhoaux(jmtd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,sqrt
C     ..
C     .. Data statements ..
      DATA czero/ (0.000,0.000)/
c
c     set components of gradient in terms of Ylms
c
      pi = pimach()
      s34 = sqrt(3.0/ (4.0*pi))
      s38 = sqrt(3.0/ (8.0*pi))
      grd1(1,0) = czero
      grd1(2,0) = czero
      grd1(3,0) = cmplx(s34,0.0)
      grd1(1,-1) = cmplx(s38,0.0)
      grd1(2,-1) = cmplx(0.0,-s38)
      grd1(3,-1) = czero
      grd1(1,1) = cmplx(-s38,0.0)
      grd1(2,1) = cmplx(0.0,-s38)
      grd1(3,1) = czero
c
      WRITE  (6,*)
      WRITE (16,*)
      spin: DO jsp = 1,jspins
         na = 1
         types: DO n = 1,ntype
            IF (l_geo(n)) THEN

            DO i = 1,3
               forc_a3(i) = czero
            END DO
c
            nd = ntypsy(na)
c
            lat_har: DO lh = 1,nlh(nd)
               lindex = llh(lh,nd)
               IF (lindex.GT.1) EXIT
               DO i = 1,3
                  gv(i) = czero
               END DO
c    sum over all m for particular symm. harmonic
               DO mem = 1,nmem(lh,nd)
                  mindex = mlh(mem,lh,nd)
                  DO i = 1,3
                     gv(i) = gv(i) + clnu(mem,lh,nd)*grd1(i,mindex)
                  END DO
               END DO
               DO ir = 1,jri(n)
                  rhoaux(ir) = rho(ir,lh,n,jsp)/ (rmsh(ir,n)**2)*
     +                         (1.0- (rmsh(ir,n)/rmt(n))**3)
               END DO
               CALL intgr3(rhoaux,rmsh(1,n),dx(n),jri(n),w)
               a3_1 = 4.0*pi/3.0*w
c+Gu
               a3_2 = vr(jri(n),lh,n,jsp)/(jspins*rmt(n))
c-Gu
c               WRITE (16,FMT=8000) a3_1,a3_2,zatom(n),rmt(n)
c 8000          FORMAT (' a3_1,a3_2,zatom(n),rmt(n)',4e12.4)
               DO i = 1,3
                  forc_a3(i) = forc_a3(i) + (a3_1+a3_2)*gv(i)*zatom(n)
               END DO
            ENDDO lat_har

            DO i = 1,3
               force(i,n,jsp) = force(i,n,jsp) + real(forc_a3(i))
            END DO
c
c     write result
            WRITE (6,FMT=8010) n
            WRITE (16,FMT=8010) n
            WRITE (6,FMT=8020) (forc_a3(i),i=1,3)
            WRITE (16,FMT=8020) (forc_a3(i),i=1,3)
 8010       FORMAT (' FORCES: EQUATION A3 FOR ATOM TYPE',i4)
 8020       FORMAT (' FX_A3=',2f10.6,' FY_A3=',2f10.6,' FZ_A3=',2f10.6)

            ENDIF               ! l_geo(n) = .true.
            na = na + neq(n)
         ENDDO types

      ENDDO spin

      END SUBROUTINE force_a3
      END MODULE m_force_a3
