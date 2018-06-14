      MODULE m_outcdn
c     ********************************************************
c     calculates the charge density at given point p(i=1,3)
c     ********************************************************
      CONTAINS
      SUBROUTINE outcdn(
     >                  p,n,na,iv,iflag,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                  n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd,jmtd,
     >                  natd,ntypd,nmzd,nop,nop2,mrot,tau,symor,invtab,
     >                  nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,mlh,llh,
     >                  clnu,nmem,nlh,lmax,rmsh,jri,pos,ngopr,ntypsy,
     >                  amat,bmat,qpw,rhtxy,rho,rht,odi,ods,
     <                  xdnout)
c
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_angle
      USE m_cotra, ONLY : cotra0,cotra1
      USE m_starf, ONLY : starf2,starf3
      USE m_ylm
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd
      INTEGER, INTENT (IN) :: lmaxd,jmtd,ntypd,natd,nmzd
      INTEGER, INTENT (IN) :: iflag,jsp,n,na,iv
      INTEGER, INTENT (IN) :: nq3,nvac,nmz,nmzxy,nq2,nop,nop2
      LOGICAL, INTENT (IN) :: plpot,symor,invs
      REAL,    INTENT (IN) :: z1,delz
      REAL,    INTENT (OUT) :: xdnout
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: qpw(n3d,jspd)
      COMPLEX, INTENT (IN) :: rhtxy(nmzxyd,odi%n2d-1,2,jspd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ngopr(natd),ntypsy(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),lmax(ntypd)
      INTEGER, INTENT (IN) :: kv2(2,n2d),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),pos(3,natd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),tau(3,nop)
      REAL,    INTENT (INOUT) :: p(3)
C     ..
C     .. Local Scalars ..
      REAL delta,s,sx,xd1,xd2,xx1,xx2,rrr,phi
      INTEGER i,j,jp3,jr,k,lh,mem,nd,nopa,ivac,ll1,lm,m,gzi
      COMPLEX ci
C     ..
C     .. Local Arrays ..
      COMPLEX sf2(n2d),sf3(n3d),ylm((lmaxd+1)**2)
      REAL rcc(3),x(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,real,sqrt
C     ..
      ci = cmplx(0.,1.)
      ivac=iv
      IF (iflag.EQ.0) GO TO 20
      IF (iflag.EQ.1) GO TO 40
c     ---> interstitial part
      CALL cotra1(p(1),rcc,bmat)
      CALL starf3(
     >            nop,nq3,symor,kv3,mrot,tau,rcc,invtab,
     <            sf3)
c
      xdnout = 0.0
      DO 10 k = 1,nq3
         xdnout = xdnout + real(qpw(k,jsp)*sf3(k))*nstr(k)
   10 CONTINUE
      RETURN
c     ---> vacuum part
   20 CONTINUE
      xdnout = 0.
c-odim
      IF (odi%d1) THEN
         rrr = sqrt( p(1)**2 + p(2)**2 )
         phi = angle(p(1),p(2))
         jp3 = (rrr-z1)/delz
         delta = (rrr-z1)/delz - jp3
*     we count 0 as point 1
         jp3 = jp3 + 1
         IF (jp3.LT.nmz) THEN
            xdnout = rht(jp3,ivac,jsp) + delta*
     +           (rht(jp3+1,ivac,jsp)-rht(jp3,ivac,jsp))
            IF (jp3.LT.nmzxy) THEN
               xx1 = 0.
               xx2 = 0.
               DO 30 k = 2,odi%nq2
                  m = odi%kv(2,k)
                  gzi = odi%kv(1,k)
                  xx1 = xx1 + real(rhtxy(jp3,k-1,ivac,jsp)*
     *                 exp(ci*m*phi)*exp(ci*gzi*bmat(3,3)*p(3)))*
     *                 odi%nst2(k)
                  xx2 = xx2 + real(rhtxy(jp3+1,k-1,ivac,jsp)*
     *                 exp(ci*m*phi)*exp(ci*gzi*bmat(3,3)*p(3)))*
     *                 odi%nst2(k)
 30            CONTINUE
               xdnout = xdnout + xx1 + delta* (xx2-xx1)
            END IF
         ELSE
            xdnout = 0.0
         END IF

      ELSE
c+odim      
         IF (p(3).LT.0.0) THEN
            ivac = nvac
            IF (invs) THEN
               p(1) = -p(1)
               p(2) = -p(2)
            END IF
            p(3) = abs(p(3))
         END IF
*   below: delete this - temporary
*      if (invs) then
*      p(1)=-p(1)
*      p(2)=-p(2)
*      endif
*
         CALL cotra1(p,rcc,bmat)
         CALL starf2(
     >            nop2,nq2,kv2,mrot,symor,tau,rcc,invtab,
     <            sf2)
c
         jp3 = (p(3)-z1)/delz
         delta = (p(3)-z1)/delz - jp3
*     we count 0 as point 1
         jp3 = jp3 + 1
         IF (jp3.LT.nmz) THEN
             xdnout = rht(jp3,ivac,jsp) + delta*
     +               (rht(jp3+1,ivac,jsp)-rht(jp3,ivac,jsp))
            IF (jp3.LT.nmzxy) THEN
               xx1 = 0.
               xx2 = 0.
              DO 35 k = 2,nq2
               xx1 = xx1 + real(rhtxy(jp3,k-1,ivac,jsp)*sf2(k))*nstr2(k)
               xx2 = xx2 + real(rhtxy(jp3+1,k-1,ivac,jsp)*sf2(k))*
     +               nstr2(k)
   35         CONTINUE
              xdnout = xdnout + xx1 + delta* (xx2-xx1)
            END IF
         ELSE
            xdnout = 0.0
         END IF
c----> vacuum finishes
      ENDIF

      RETURN
c     ----> m.t. part
   40 CONTINUE
      
      nd = ntypsy(na)
      nopa = ngopr(na)
      IF (odi%d1) nopa = ods%ngopr(na)
      sx = 0.0
      DO 50 i = 1,3
         x(i) = p(i) - pos(i,na)
         sx = sx + x(i)*x(i)
   50 CONTINUE
      sx = sqrt(sx)
      IF (nopa.NE.1) THEN
c... switch to internal units
         CALL cotra1(x,rcc,bmat)
c... rotate into representative
         DO 70 i = 1,3
            p(i) = 0.
            DO 60 j = 1,3
              IF (.NOT.odi%d1) THEN
               p(i) = p(i) + mrot(i,j,nopa)*rcc(j)
              ELSE
               p(i) = p(i) + ods%mrot(i,j,nopa)*rcc(j)
              END IF
   60       CONTINUE
   70    CONTINUE
c... switch back to cartesian units
         CALL cotra0(p,x,amat)
      END IF
      DO 80 j = jri(n),2,-1
         IF (sx.GE.rmsh(j,n)) GO TO 90
   80 CONTINUE
   90 jr = j
      CALL ylm4(
     >          lmax(n),x,
     <          ylm)
      xd1 = 0.0
      xd2 = 0.0
      DO 110 lh = 0, nlh(nd)
         ll1 = llh(lh,nd) * ( llh(lh,nd) + 1 ) + 1
         s = 0.0
         DO mem = 1,nmem(lh,nd)
           lm = ll1 + mlh(mem,lh,nd)
           s = s + real( clnu(mem,lh,nd)*ylm(lm) )
         ENDDO
         IF (plpot) THEN
            xd1 = xd1 + rho(jr,lh,n,jsp)*s
         ELSE
            xd1 = xd1 + rho(jr,lh,n,jsp)*s/ (rmsh(jr,n)*rmsh(jr,n))
         END IF
         IF (jr.EQ.jri(n)) GO TO 110
         IF (plpot) THEN
            xd2 = xd2 + rho(jr+1,lh,n,jsp)*s
         ELSE
            xd2 = xd2 + rho(jr+1,lh,n,jsp)*s/
     +            (rmsh(jr+1,n)*rmsh(jr+1,n))
         END IF
  110 CONTINUE
      IF (jr.EQ.jri(n)) THEN
         xdnout = xd1
      ELSE
         xdnout = xd1 + (xd2-xd1) *
     +                  (sx-rmsh(jr,n)) / (rmsh(jr+1,n)-rmsh(jr,n))
      END IF
 8000 FORMAT (2f10.6)
c
      RETURN
      END SUBROUTINE outcdn
      END MODULE m_outcdn
