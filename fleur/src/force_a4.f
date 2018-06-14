      MODULE m_force_a4
      CONTAINS
      SUBROUTINE force_a4(jmtd,memd,nlhd,ntypd,jspd,ntypsd,natd,
     >                    jri,mlh,nlh,nmem,ntype,jspins,ntypsy,
     >                    dx,rmsh,clnu,llh,vr,neq,l_geo,
     X                    force)
c
c ************************************************************
c rhocore force contribution a la Rici et al.
c
c ************************************************************
c
      USE m_intgr, ONLY : intgr0,intgr3
      USE m_constants, ONLY : pimach
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,memd,nlhd,ntypd,jspd,ntypsd
      INTEGER, INTENT (IN) :: ntype,jspins,natd
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nlh(ntypsd),mlh(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd,ntypd,jspd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
C     ..
C     .. Local Scalars ..
      COMPLEX czero
      REAL a4_1,a4_2,qcore,s13,s23,w,xi,sfp
      INTEGER i,ir,jsp,lh,lindex,mem,mindex,n,nd,na
C     ..
C     .. Local Arrays ..
      COMPLEX forc_a4(3),gv(3),ycomp1(3,-1:1)
      REAL rhoaux(jmtd),rhoc(jmtd)
C     ..
C     .. External Functions ..
      REAL difcub
      EXTERNAL difcub
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,sqrt
C     ..
C     .. Data statements ..
      DATA czero/ (0.000,0.000)/
C     ..
c
c     set ylm related components
c
      sfp = 2 * sqrt(pimach())
      s13 = sqrt(1.0/3.0)
      s23 = sqrt(2.0/3.0)
      ycomp1(1,0) = czero
      ycomp1(2,0) = czero
      ycomp1(3,0) = cmplx(2.0*s13,0.0)
      ycomp1(1,-1) = cmplx(s23,0.0)
      ycomp1(2,-1) = cmplx(0.0,-s23)
      ycomp1(3,-1) = czero
      ycomp1(1,1) = cmplx(-s23,0.0)
      ycomp1(2,1) = cmplx(0.0,-s23)
      ycomp1(3,1) = czero
c
      OPEN (17,file='cdnc',form='unformatted',status='old')
      DO 40 jsp = 1,jspins
         IF (jsp.EQ.1) THEN
            REWIND 17
         END IF
         na = 1
         DO 30 n = 1,ntype
c     --->    read in core density
            READ (17) (rhoc(ir),ir=1,jri(n))
c     --->   core density is multiplied by r**2*sqrt(4*pi)
c     --->   such that intrg0 gives correct core charge
c
c     --->    read in kinetic enrgy of the core
            READ (17)
c
            IF (l_geo(n)) THEN
            DO i = 1,3
               forc_a4(i) = czero
            END DO
c
            CALL intgr0(rhoc,rmsh(1,n),dx(n),jri(n),qcore)
c     write(16,1616) qcore
 8000       FORMAT (' FORCE_A4: core charge=',1p,e16.8)
c
c
            nd = ntypsy(na)
c
c
            DO 10 lh = 1,nlh(nd)
               lindex = llh(lh,nd)
               IF (lindex.GT.1) GO TO 20
               DO i = 1,3
                  gv(i) = czero
               END DO
c
c    sum over all m for particular symm. harmonic
               DO mem = 1,nmem(lh,nd)
                  mindex = mlh(mem,lh,nd)
                  DO i = 1,3
                     gv(i) = gv(i) + clnu(mem,lh,nd)*ycomp1(i,mindex)
                  END DO
               END DO
c
c
c     construct integrand rhocore*d/dr(v)*r**2
c     note: rhocore is already multiplied by r**2 and srt(4.*pi)
c     difcub performs analytic derivative of Lagrangian of 3rd order
c
               xi = rmsh(1,n)
               rhoaux(1) = difcub(rmsh(1,n),vr(1,lh,n,jsp),xi)*rhoc(1)
               DO ir = 2,jri(n) - 2
                  xi = rmsh(ir,n)
                  rhoaux(ir) = difcub(rmsh(ir-1,n),
     +                                vr(ir-1,lh,n,jsp),xi) * rhoc(ir)
               END DO
c
               ir = jri(n) - 1
               xi = rmsh(ir,n)
               rhoaux(ir) = difcub(rmsh(jri(n)-3,n),
     +                              vr(jri(n)-3,lh,n,jsp),xi)*rhoc(ir)
c
               ir = jri(n)
               xi = rmsh(ir,n)
               rhoaux(ir) = difcub(rmsh(jri(n)-3,n),
     +                              vr(jri(n)-3,lh,n,jsp),xi)*rhoc(ir)
               CALL intgr3(rhoaux,rmsh(1,n),dx(n),jri(n),w)
               a4_1 = 0.5*w/sfp
c
c     construct integrand rhocore*v*r
c     note: rhocore is already multiplied by r**2 and srt(4.*pi)
c
               DO ir = 1,jri(n)
                  rhoaux(ir) = rhoc(ir)/rmsh(ir,n)*vr(ir,lh,n,jsp)
               END DO
c
               CALL intgr3(rhoaux,rmsh(1,n),dx(n),jri(n),w)
               a4_2 = w/sfp
c
               DO i = 1,3
                  forc_a4(i) = forc_a4(i) - (a4_1+a4_2)*gv(i)
               END DO
c
c  lh loop ends
   10       CONTINUE
   20       CONTINUE
c
c
c     sum to existing forces
c
            DO i = 1,3
               force(i,n,jsp) = force(i,n,jsp) + real(forc_a4(i))
            END DO
c
c     write result
c
            WRITE (6,FMT=8010) n
            WRITE (16,FMT=8010) n
            WRITE (6,FMT=8020) (forc_a4(i),i=1,3)
            WRITE (16,FMT=8020) (forc_a4(i),i=1,3)
 8010       FORMAT (' FORCES: EQUATION A4 FOR ATOM TYPE',i4)
 8020       FORMAT (' FX_A4=',2f10.6,' FY_A4=',2f10.6,' FZ_A4=',2f10.6)
c type loop ends
            ENDIF
            na = na + neq(n)
   30    CONTINUE
         READ (17)
c spin loop ends
   40 CONTINUE
c
      CLOSE (17)
      RETURN
      END SUBROUTINE force_a4
      END MODULE m_force_a4
