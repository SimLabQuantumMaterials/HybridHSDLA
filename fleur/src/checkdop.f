      MODULE m_checkdop
      CONTAINS
      SUBROUTINE checkdop(
     >                    p,np,n,na,ivac,iflag,jsp,cdn,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    fpw,fr,fxy,fz,odi,ods)
c ************************************************************
c     subroutines checks the continuity of coulomb           *
c     potential or valence charge density                    *
c     across the m.t. and vacuum boundaries                  *
c                                     c.l.fu                 *
c     (unifies vcheck and cdnchk      g.b.                   *
c YM:  this routine doesn't really work in the vacuum in 1D case yet
c ************************************************************

      USE m_cotra, ONLY : cotra0,cotra1
      USE m_starf, ONLY : starf2,starf3
      USE m_angle
      USE m_ylm
      USE m_od_types, ONLY : od_inp,od_sym

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,n2d,n3d,ntypd
      INTEGER, INTENT (IN) :: lmaxd,jspd,natd,nmzd,nmzxyd,nspd,jmtd
      INTEGER, INTENT (IN) :: iflag,ivac,n,na,np,jsp,nop,nop2,nq2,nq3
      LOGICAL, INTENT (IN) :: cdn,symor
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: kv2(2,n2d),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),jri(ntypd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nlh(ntypsd),ngopr(natd)
      INTEGER, INTENT (IN) :: ntypsy(natd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),lmax(ntypd)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),tau(3,nop)
      REAL,    INTENT (IN) :: p(3,nspd),rmt(ntypd),pos(3,natd)
      REAL,    INTENT (IN) :: fr(jmtd,0:nlhd,ntypd,jspd),fz(nmzd,2,jspd)
      COMPLEX, INTENT (IN) :: fpw(n3d,jspd),fxy(nmzxyd,odi%n2d-1,2,jspd)
C     ..
C     .. Local Scalars ..
      REAL av,dms,rms,s,ir2,help,phi
      INTEGER i,j,k,lh,mem,nd,lm,ll1,nopa,m,gz
      COMPLEX ic
C     ..
C     .. Local Arrays ..
      COMPLEX sf2(n2d),sf3(n3d),ylm( (lmaxd+1)**2 )
      REAL rcc(3),v1(nspd),v2(nspd),x(3),ri(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL fitchk
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
      ic = cmplx(0.,1.)
      IF (iflag.LT.0) GO TO 50
c     ---> Interstitial part
      DO j = 1,np
         IF (.NOT.odi%d1) THEN
           CALL cotra1(p(:,j),rcc,bmat)
           CALL starf3(
     >               nop,nq3,symor,kv3,mrot,tau,p(1,j),invtab,
     <               sf3)
         ENDIF
c
         IF (odi%d1) THEN
           CALL cotra1(p(:,j),rcc,bmat)
           CALL starf3(
     >               nop,nq3,symor,kv3,mrot,tau,rcc,invtab,
     <               sf3)
         ENDIF
         v1(j) = 0.0
         DO k = 1,nq3
            v1(j) = v1(j) + real(fpw(k,jsp)*sf3(k))*nstr(k)
         ENDDO
      ENDDO
c     ---> vacuum part
      IF (cdn) THEN
        WRITE (6,FMT=9000) ivac
        WRITE (16,FMT=9000) ivac
      ELSE
        WRITE (6,FMT=8000) ivac
        WRITE (16,FMT=8000) ivac
      ENDIF
      DO 40 j = 1,np
         IF (.NOT.odi%d1) THEN
            CALL starf2(
     >           nop2,nq2,kv2,mrot,symor,tau,p(1,j),invtab,
     <           sf2)
            v2(j) = fz(1,ivac,jsp)
            DO 30 k = 2,nq2
               v2(j) = v2(j) + real(fxy(1,k-1,ivac,jsp)*sf2(k))*nstr2(k)
   30       CONTINUE
         ELSE
c-odim
            v2(j) = fz(1,ivac,jsp)
            phi = angle(p(1,j),p(2,j))
            DO 35 k = 2,odi%nq2
               m = odi%kv(2,k)
               gz = odi%kv(1,k)
               v2(j) = v2(j) + real(fxy(1,k-1,ivac,jsp)*
     *           exp(ic*m*phi)*exp(ic*bmat(3,3)*gz*p(3,j)))*odi%nst2(k)
 35         CONTINUE
c+odim
         END IF
         IF (odi%d1) THEN
            CALL cotra1(p(:,j),rcc,bmat)
            WRITE (6,FMT=8020)  rcc,(p(i,j),i=1,3),v1(j),v2(j)
            WRITE (16,FMT=8020) rcc,(p(i,j),i=1,3),v1(j),v2(j)
         ELSE
            CALL cotra0(p(1,j),rcc,amat)
            WRITE (6,FMT=8020) (p(i,j),i=1,3),rcc,v1(j),v2(j)
            WRITE (16,FMT=8020) (p(i,j),i=1,3),rcc,v1(j),v2(j)
         ENDIF
   40 CONTINUE
      CALL fitchk(v1,v2,np,av,rms,dms)
      WRITE (6,FMT=8030) av,rms,dms
      WRITE (16,FMT=8030) av,rms,dms
      RETURN
c      ----> interstitial part
   50 CONTINUE
      DO j = 1,np
         CALL cotra1(p(1,j),rcc,bmat)
         CALL starf3(
     >               nop,nq3,symor,kv3,mrot,tau,rcc,invtab,
     <               sf3)
c
         v1(j) = 0.0
         DO k = 1,nq3
            v1(j) = v1(j) + real(fpw(k,jsp)*sf3(k))*nstr(k)
         ENDDO
      ENDDO
c     ----> m.t. part
      IF (cdn) THEN
        WRITE (6,FMT=9010) n
        WRITE (16,FMT=9010) n
      ELSE
        WRITE (6,FMT=8010) n
        WRITE (16,FMT=8010) n
      ENDIF
      ir2 = 1.0
      IF (cdn) ir2 = 1.0 / ( rmt(n)*rmt(n) )
      nd = ntypsy(na)
      nopa = ngopr(na)
      IF (odi%d1) THEN
         nopa = ods%ngopr(na)
         nopa = ods%invtab(nopa)
      END IF
      DO 110 j = 1,np
         DO i = 1,3
            x(i) = p(i,j) - pos(i,na)
         ENDDO
! new
         IF (nopa.NE.1) THEN
            CALL cotra1(x,rcc,bmat)  ! switch to internal units
            DO i = 1,3               ! rotate into representative
               ri(i) = 0.
               DO k = 1,3
                  IF (odi%d1) THEN
                     ri(i) = ri(i) + ods%mrot(i,k,nopa)*rcc(k)
                  ELSE
                     ri(i) = ri(i) + mrot(i,k,nopa)*rcc(k)
                  END IF
               ENDDO
             ENDDO
            CALL cotra0(ri,x,amat)    !switch back to cartesian units
         END IF
! new
         CALL ylm4(
     >             lmax(n),x,
     <             ylm)
         help = 0.0
         DO lh = 0,nlh(nd)
            s = 0.0
            ll1 = llh(lh,nd) * ( llh(lh,nd) + 1 ) + 1
            DO mem = 1,nmem(lh,nd)
               lm = ll1 + mlh(mem,lh,nd)
               s = s + real( clnu(mem,lh,nd)* ylm(lm) )
            ENDDO
            help = help + fr(jri(n),lh,n,jsp) * s
         ENDDO
         v2(j) = help * ir2 
         IF (j.LE.8) THEN
            CALL cotra1(p(1,j),rcc,bmat)
            WRITE (6,FMT=8020) rcc, (p(i,j),i=1,3),v1(j),v2(j)
            WRITE (16,FMT=8020) rcc, (p(i,j),i=1,3),v1(j),v2(j)
         END IF
  110 CONTINUE
      CALL fitchk(v1,v2,np,av,rms,dms)
      WRITE (6,FMT=8030) av,rms,dms
      WRITE (16,FMT=8030) av,rms,dms
 8000 FORMAT (/,'    int.-vac. boundary (potential): ivac=',i2,/,t10,
     +       'int-coord',t36,'cart-coord',t57,' inter. ',t69,' vacuum ')
 8010 FORMAT (/,'    int.-m.t. boundary (potential): atom type=',i2,/,
     +       t10,'int-coord',t36,'cart-coord',t57,' inter. ',t69,
     +       '  m. t. ')
 8020 FORMAT (1x,2 (3f8.3,2x),2f12.6)
 8030 FORMAT (/,10x,'average value = ',f10.6,/,t11,'rms,dmx=',2f7.3,
     +       ' per cent')
 9000 FORMAT (/,'    int.-vac. boundary (density): ivac=',i2,/,t10,
     +       'int-coord',t36,'cart-coord',t57,' inter. ',t69,' vacuum ')
 9010 FORMAT (/,'    int.-m.t. boundary (density): atom type=',i2,/,t10,
     +       'int-coord',t36,'cart-coord',t57,' inter. ',t69,'  m. t. ')
      RETURN
      END SUBROUTINE checkdop
      END MODULE m_checkdop
