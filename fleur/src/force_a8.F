      MODULE m_forcea8
c ************************************************************
c Pulay 1st term  force contribution a la Rici et al., eq. A8
c
c ************************************************************
      CONTAINS
      SUBROUTINE force_a8(
     >                    jmtd,memd,nlhd,ntypd,jspd,ntypsd,lmaxd,l_geo,
     >                    jri,mlh,nlh,nmem,ntype,jsp,ntypsy,natd,
     >                    dx,rmsh,clnu,llh,vr,rho,neq,
     >                    f_a12,f_a21,f_b4,f_b8,
     X                    force)
c
      USE m_intgr, ONLY : intgr3
      USE m_constants, ONLY : pimach
      USE m_gaunt, ONLY :gaunt1
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,memd,nlhd,ntypd,jspd,ntypsd,lmaxd
      INTEGER, INTENT (IN) :: ntype,jsp,natd
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nlh(ntypsd),mlh(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: f_a12(3,ntypd),f_a21(3,ntypd)
      COMPLEX, INTENT (IN) :: f_b4(3,ntypd),f_b8(3,ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd,ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
C     ..
C     .. Local Scalars ..
      COMPLEX aaa,bbb,ccc,czero,ddd,eee,fff,ci
      REAL a8_1,a8_2,qval,rr,truerho,truev,xi,sfp,pi
      INTEGER i,ir,j,l,l1,l2,lh1,lh2,m,m1,m2,mem1,mem2,n,nd,na
C     ..
C     .. Local Arrays ..
      COMPLEX forc_a8(3),gv(3),f_sum(3)
      REAL rhoaux(jmtd),rhodif(jmtd)
C     ..
C     .. External Functions ..
      REAL difcub
      EXTERNAL difcub
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min,real,sign,sqrt,cmplx
C     ..
C     .. Statement Functions ..
      REAL alpha,beta,delta,epslon,gamma,phi
      INTEGER krondel
C     ..
C     .. Data statements ..
      DATA czero/ (0.000,0.000)/
c
c  inline functions:
c
c Kronecker delta for arguments >=0 AND <0
c
c
      krondel(i,j) = min(abs(i)+1,abs(j)+1)/max(abs(i)+1,abs(j)+1)*
     +               (1+sign(1,i)*sign(1,j))/2
      alpha(l,m) = (l+1)*0.5e0*sqrt(real((l-m)* (l-m-1))/
     +             real((2*l-1)* (2*l+1)))
      beta(l,m) = l*0.5e0*sqrt(real((l+m+2)* (l+m+1))/
     +            real((2*l+1)* (2*l+3)))
      gamma(l,m) = (l+1)*0.5e0*sqrt(real((l+m)* (l+m-1))/
     +             real((2*l-1)* (2*l+1)))
      delta(l,m) = l*0.5e0*sqrt(real((l-m+2)* (l-m+1))/
     +             real((2*l+1)* (2*l+3)))
      epslon(l,m) = (l+1)*sqrt(real((l-m)* (l+m))/
     +               real((2*l-1)* (2*l+1)))
      phi(l,m) = l*sqrt(real((l-m+1)* (l+m+1))/real((2*l+1)* (2*l+3)))
C     ..
      ci = cmplx(0.0,1.0)
      pi = pimach()
      sfp = sqrt( 4*pi )
      WRITE  (6,*)
      WRITE (16,*)
c
      na = 1
      DO 50 n = 1,ntype
         IF (l_geo(n)) THEN
         DO i = 1,3
            forc_a8(i) = czero
         END DO
c
c
         nd = ntypsy(na)
c
         CALL intgr3(rho(1,0,n,jsp),rmsh(1,n),dx(n),jri(n),qval)
c
c     check if l=0 density is correct;
c     note that in general also all l>0
c     components of the density have been multiplied by r**2
c     (see for example subr. checkdop which constructs true MT-density at Rmt)
c     factor sqrt(4pi) comes from Y_00 * \int d\Omega = 1/sqrt(4pi) * 4pi
c     write(16,1616) qval*sfp
 8000    FORMAT (' FORCE_A8: valence charge=',1p,e16.8)
c
c    PART I of FORCE_A8
c
         DO 20 lh1 = 0,nlh(nd)
            l1 = llh(lh1,nd)
            DO 10 lh2 = 0,nlh(nd)
               l2 = llh(lh2,nd)
c
               DO i = 1,3
                  gv(i) = czero
               END DO
c
c    sum over all m for particular symm. harmonic
               DO mem1 = 1,nmem(lh1,nd)
                  m1 = mlh(mem1,lh1,nd)
                  DO mem2 = 1,nmem(lh2,nd)
                     m2 = mlh(mem2,lh2,nd)
                     gv(1) = gv(1) + sqrt(2.e0*pi/3.e0)*
     +                       clnu(mem1,lh1,nd)*clnu(mem2,lh2,nd)*
     +                       (gaunt1(1,l1,l2,-1,m1,m2,lmaxd)-
     +                       gaunt1(1,l1,l2,1,m1,m2,lmaxd))
                     gv(2) = gv(2) - ci*sqrt(2.e0*pi/3.e0)*
     +                       clnu(mem1,lh1,nd)*clnu(mem2,lh2,nd)*
     +                       (gaunt1(1,l1,l2,-1,m1,m2,lmaxd)+
     +                       gaunt1(1,l1,l2,1,m1,m2,lmaxd))
                     gv(3) = gv(3) + sqrt(4.e0*pi/3.e0)*
     +                       clnu(mem1,lh1,nd)*clnu(mem2,lh2,nd)*
     +                       gaunt1(1,l1,l2,0,m1,m2,lmaxd)
                  END DO
               END DO
c
c     note that in general also all l>0
c     components of the density have been multiplied by r**2
c     here we need the true radial denisity for performing the derivative
c     therefore we divide by r**2
               DO ir = 1,jri(n)
                  rhoaux(ir) = rho(ir,lh2,n,jsp)/ (rmsh(ir,n)**2)
               END DO
c
c
C NOTE: we should have here: vr  = vtrue
c     difcub performs analytic derivative of Lagrangian of 3rd order
c
               xi = rmsh(1,n)
               rr = xi*xi
               rhodif(1) = difcub(rmsh(1,n),rhoaux(1),xi)*
     +                     vr(1,lh1,n)*rr
               DO ir = 2,jri(n) - 2
                  xi = rmsh(ir,n)
                  rr = xi*xi
                  rhodif(ir) = difcub(rmsh(ir-1,n),rhoaux(ir-1),xi)*
     +                         vr(ir,lh1,n)*rr
               END DO
c
               xi = rmsh(jri(n)-1,n)
               rr = xi*xi
               rhodif(jri(n)-1) = difcub(rmsh(jri(n)-3,n),
     +                            rhoaux(jri(n)-3),xi)*
     +                            vr(jri(n)-1,lh1,n)*rr
c
               xi = rmsh(jri(n),n)
               rr = xi*xi
               rhodif(jri(n)) = difcub(rmsh(jri(n)-3,n),
     +                          rhoaux(jri(n)-3),xi)*
     +                          vr(jri(n),lh1,n)*rr
C NOTE: vr(l=0) is EXPLICITELY multiplied by r/sqrt(4pi) to be the TRUE
C r*V which is needed for the radial Schr. equ.
C here, we need the l=0 component, v_l=0(r) which will be multiplied by
C Y_00 in the lm expansion; therefore we MUST recorrect vr(l=0) by
C the incerse factor sqrt(4pi)/r. We do the correction for the product
C array rhodif
               IF (lh1.EQ.0) THEN
                  DO ir = 1,jri(n)
                     rhodif(ir) = rhodif(ir)/rmsh(ir,n)*sfp
                  END DO
               END IF
               CALL intgr3(rhodif,rmsh(1,n),dx(n),jri(n),a8_1)
c
c
               DO i = 1,3
                  forc_a8(i) = forc_a8(i) + a8_1*gv(i)
               END DO
c
c  lh1,lh2 loops end
   10       CONTINUE
   20    CONTINUE
c     END OF PART I
c
c    PART II of FORCE_A8
c
         DO 40 lh1 = 0,nlh(nd)
            l1 = llh(lh1,nd)
            DO 30 lh2 = 0,nlh(nd)
               l2 = llh(lh2,nd)
c
               DO ir = 1,jri(n)
                  truev = vr(ir,lh1,n)
                  truerho = rho(ir,lh2,n,jsp)/ (rmsh(ir,n)**2)
                  rhoaux(ir) = truev*truerho*rmsh(ir,n)
               END DO
               IF (lh1.EQ.0) THEN
                  DO ir = 1,jri(n)
                     rhoaux(ir) = rhoaux(ir)/rmsh(ir,n)*sfp
                  END DO
               END IF
               CALL intgr3(rhoaux,rmsh(1,n),dx(n),jri(n),a8_2)
c
               DO i = 1,3
                  gv(i) = (0.0,0.0)
               END DO
c
c    sum over all m for particular sym. harmonic
               DO mem1 = 1,nmem(lh1,nd)
                  m1 = mlh(mem1,lh1,nd)
                  DO mem2 = 1,nmem(lh2,nd)
                     m2 = mlh(mem2,lh2,nd)
C  NOTE: delta(-m,m')(-1)^m was applied because we have integrand
C        Y nabla Y instead of Y* nabla Y.
C        Then we know that  Y*(lm) = (-1)^m Y(l,-m).
C                      and  Y(lm)  = (-1)^m Y*(l,-m).
C        Therefore  (-1)^m delta (-m,m') appears
C
                     aaa = alpha(l2,m2)*clnu(mem1,lh1,nd)*
     +                     clnu(mem2,lh2,nd)* (-1)**m1*krondel(l1,l2-1)*
     +                     krondel(-m1,m2+1)
                     bbb = beta(l2,m2)*clnu(mem1,lh1,nd)*
     +                     clnu(mem2,lh2,nd)* (-1)**m1*krondel(l1,l2+1)*
     +                     krondel(-m1,m2+1)
                     ccc = gamma(l2,m2)*clnu(mem1,lh1,nd)*
     +                     clnu(mem2,lh2,nd)* (-1)**m1*krondel(l1,l2-1)*
     +                     krondel(-m1,m2-1)
                     ddd = delta(l2,m2)*clnu(mem1,lh1,nd)*
     +                     clnu(mem2,lh2,nd)* (-1)**m1*krondel(l1,l2+1)*
     +                     krondel(-m1,m2-1)
                     eee = epslon(l2,m2)*clnu(mem1,lh1,nd)*
     +                     clnu(mem2,lh2,nd)* (-1)**m1*krondel(l1,l2-1)*
     +                     krondel(-m1,m2)
                     fff = phi(l2,m2)*clnu(mem1,lh1,nd)*
     +                     clnu(mem2,lh2,nd)* (-1)**m1*krondel(l1,l2+1)*
     +                     krondel(-m1,m2)
c
                     gv(1) = gv(1) + aaa + bbb - ccc - ddd
                     gv(2) = gv(2) - ci* (aaa+bbb+ccc+ddd)
                     gv(3) = gv(3) + eee - fff
c
c  end of summation m1,m2
                  END DO
               END DO
c
               DO i = 1,3
                  forc_a8(i) = forc_a8(i) + a8_2*gv(i)
               END DO
c
c  lh1,lh2 loops end
   30       CONTINUE
   40    CONTINUE
c
c     sum to existing forces
c
         DO i = 1,3
            force(i,n,jsp) = force(i,n,jsp) + real(forc_a8(i))
         END DO
c
c     write result
c
         WRITE (6,FMT=8010) n
         WRITE (16,FMT=8010) n
         WRITE (6,FMT=8020) (forc_a8(i),i=1,3)
         WRITE (16,FMT=8020) (forc_a8(i),i=1,3)
 8010    FORMAT (' FORCES: EQUATION A8 FOR ATOM TYPE',i4)
 8020    FORMAT (' FX_A8=',2f10.6,' FY_A8=',2f10.6,' FZ_A8=',2f10.6)

         ENDIF
         na = na + neq(n)
   50 CONTINUE
c
c
c     write now also result of a12 & a21 ( also b4 and b8 )
c     and b4
#ifndef CPP_APW
      WRITE  (6,*)
      WRITE (16,*)
      DO n=1,ntype
         IF (l_geo(n)) THEN
            WRITE  (6,FMT=8030) n
            WRITE (16,FMT=8030) n
            WRITE  (6,FMT=8040) (f_a12(i,n),i=1,3)
            WRITE (16,FMT=8040) (f_a12(i,n),i=1,3)
         ENDIF
 8030    FORMAT (' FORCES: EQUATION A12 FOR ATOM TYPE',i4)
 8040    FORMAT (' FX_A12=',2f10.6,' FY_A12=',2f10.6,' FZ_A12=',2f10.6)
      ENDDO
#else
      WRITE  (6,*)
      WRITE (16,*)
      DO n=1,ntype
         IF (l_geo(n)) THEN
            WRITE  (6,FMT=8070) n
            WRITE (16,FMT=8070) n
            WRITE  (6,FMT=8080) (f_b4(i,n),i=1,3)
            WRITE (16,FMT=8080) (f_b4(i,n),i=1,3)
         ENDIF
 8070    FORMAT (' FORCES: EQUATION B4 FOR ATOM TYPE',i4)
 8080    FORMAT (' FX_B4=',2f10.6,' FY_B4=',2f10.6,' FZ_B4=',2f10.6)
      ENDDO
      WRITE  (6,*)
      WRITE (16,*)
      DO n=1,ntype
         IF (l_geo(n)) THEN
            WRITE  (6,FMT=8090) n
            WRITE (16,FMT=8090) n
            WRITE  (6,FMT=8100) (f_b8(i,n),i=1,3)
            WRITE (16,FMT=8100) (f_b8(i,n),i=1,3)
         ENDIF
 8090    FORMAT (' FORCES: EQUATION B8 FOR ATOM TYPE',i4)
 8100    FORMAT (' FX_B8=',2f10.6,' FY_B8=',2f10.6,' FZ_B8=',2f10.6)
      ENDDO
#endif
      WRITE  (6,*)
      WRITE (16,*)
      DO n=1,ntype
         IF (l_geo(n)) THEN
            WRITE  (6,FMT=8050) n
            WRITE (16,FMT=8050) n
            WRITE  (6,FMT=8060) (f_a21(i,n),i=1,3)
            WRITE (16,FMT=8060) (f_a21(i,n),i=1,3)
         ENDIF
 8050    FORMAT (' FORCES: EQUATION A21 FOR ATOM TYPE',i4)
 8060    FORMAT (' FX_A21=',2f10.6,' FY_A21=',2f10.6,' FZ_A21=',2f10.6)
      ENDDO

      END SUBROUTINE force_a8
      END MODULE m_forcea8
