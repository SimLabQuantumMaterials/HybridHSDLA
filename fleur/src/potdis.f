      MODULE m_potdis
      CONTAINS
      SUBROUTINE potdis(
     >                  jspd,n3d,n2d,nmzd,nmzxyd,ntypd,nlhd,ntypsd,natd,
     >                  jspins,nq3,nq2,nmz,nmzxy,ntype,nlh,ntypsy,
     >                  jmtd,jri,dx,rmsh,neq,omtil,area,vol,
     >                  film,nvac,zrfs,invs,invs2,delz,k1d,k2d,k3d,
     >                  nk1,nk2,nk3,mx3,mx2,mx1,ig2,rgphs,ig,ustep)
c
c     *****************************************************
c     calculates the root-mean-square distance between input
c     and output potentials (averaged over unit cell volume)
c                                                     c.l.fu
c     *****************************************************
c
      USE m_intgr, ONLY : intgr3, intgz0
      USE m_constants, ONLY : pimach
      USE m_loddop
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,n3d,n2d,nmzd,nmzxyd,ntypd,nlhd,ntypsd
      INTEGER, INTENT (IN) :: jspins,nq3,nq2,nmz,nmzxy,ntype,jmtd,nvac
      INTEGER, INTENT (IN) :: nk1,nk2,nk3,mx3,mx2,mx1,k1d,k2d,k3d,natd
      LOGICAL, INTENT (IN) :: film,invs,invs2,zrfs
      REAL,    INTENT (IN) :: delz,omtil,area,vol
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),nlh(ntypsd),ntypsy(natd)
      INTEGER, INTENT (IN) :: neq(ntypd),ig2(n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),dx(ntypd)
      COMPLEX, INTENT (IN) :: ustep(n3d)
C     ..
C     .. Local Scalars ..
      REAL fact,facv,phase,rhs,sumis,sumz,fpi
      INTEGER i,i1,i2,i3,id2,id3,io,ip,iter,ivac,j,k1,k2,k3,lh,n,
     +        nk12,npz,nt,num,na
      LOGICAL tail
      CHARACTER*8 dop,iop
C     ..
C     .. Local Arrays ..
      COMPLEX rhpw(n3d,2,2),rhv1(nmzxyd,n2d-1,2,2,2)
      REAL rhsp(jmtd,0:nlhd,ntypd,2,2),rhv0(nmzd,2,2,2)
      REAL dis(2),disz(nmzd),rh(jmtd)
      REAL af3((2*k1d+1)*(2*k2d+1)*(2*k3d+1),jspd),
     +     bf3((2*k1d+1)*(2*k2d+1)*(2*k3d+1),jspd),
     +     bf2((2*k1d+1)*(2*k2d+1),jspd),
     +     af2((2*k1d+1)*(2*k2d+1),jspd)
      CHARACTER*8 name(10)
C     ..
C     .. External Subroutines ..
      EXTERNAL cdndif,cfft
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aimag,cmplx,conjg,real,sqrt
c
      fpi = 4 * pimach()
C 
      DO 10 num = 1,2
         dis(num) = 0.
   10 CONTINUE
      tail = .true.
      npz = nmz + 1
      nt = nk1*nk2*nk3
      nk12 = nk1*nk2
      fact = omtil/real(nt)
      facv = 1.0
c     ---> reload potentials
      REWIND 9
      DO 20 io = 1,2
         CALL loddop(
     >               jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >               jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >               nlh,jri,ntypsd,ntypsy,9,natd,neq,
     <               iop,dop,iter,rhsp(1,0,1,1,io),rhpw(1,1,io),
     <               rhv0(1,1,1,io),rhv1(1,1,1,1,io),name)
   20 CONTINUE
      CLOSE (9)
      IF (jspins.EQ.1) THEN
c       ---> total potential difference
         CALL cdndif(1,1,1,1,1,2,1.0,-1.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
      ELSE
c       ---> total input potential
         CALL cdndif(1,1,2,1,1,1,1.0,1.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
c       ---> total output potential
         CALL cdndif(1,1,2,2,2,2,1.0,1.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
c       ---> input 'spin' potential
         CALL cdndif(2,1,2,1,1,1,1.0,-2.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
c       ---> output 'spin' potential
         CALL cdndif(2,1,2,2,2,2,1.0,-2.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
c       ---> potential difference
         CALL cdndif(1,1,1,1,1,2,1.0,-1.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
c       ---> spin potential difference
         CALL cdndif(2,2,2,1,1,2,1.0,-1.0,film,natd,neq,
     >               n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >               nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X               rhsp,rhpw,rhv0,rhv1)
      END IF
      DO 230 num = 1,jspins
c     ----> m.t. part
         na = 1
         DO 60 n = 1,ntype
            DO 30 j = 1,jri(n)
               rh(j) = rhsp(j,0,n,num,1)*rhsp(j,0,n,num,1)*fpi
   30       CONTINUE
            DO 50 lh = 1,nlh(ntypsy(na))
               DO 40 j = 1,jri(n)
                  rh(j) = rh(j) + rhsp(j,lh,n,num,1)*rhsp(j,lh,n,num,1)*
     +                    rmsh(j,n)*rmsh(j,n)
   40          CONTINUE
   50       CONTINUE
            CALL intgr3(rh,rmsh(1,n),dx(n),jri(n),rhs)
            dis(num) = dis(num) + rhs*neq(n)
            na = na + neq(n)
   60    CONTINUE
c     ----> interstitial part
c         ---> create density in the real space
         i = 0
         DO 90 i3 = 0,nk3 - 1
            k3 = i3
            IF (k3.GT.mx3) k3 = k3 - nk3
            DO 80 i2 = 0,nk2 - 1
               k2 = i2
               IF (k2.GT.mx2) k2 = k2 - nk2
               DO 70 i1 = 0,nk1 - 1
                  k1 = i1
                  IF (k1.GT.mx1) k1 = k1 - nk1
                  i = i + 1
                  id3 = ig(k1,k2,k3)
                  phase = rgphs(k1,k2,k3)
                  IF (id3.EQ.0) THEN
                     af3(i,1) = 0.
                     bf3(i,1) = 0.
                  ELSE
                     af3(i,1) = real(rhpw(id3,num,1))*phase
                     bf3(i,1) = aimag(rhpw(id3,num,1))*phase
                  END IF
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         CALL cfft(af3(1,1),bf3(1,1),nt,nk1,nk1,1)
         CALL cfft(af3(1,1),bf3(1,1),nt,nk2,nk12,1)
         CALL cfft(af3(1,1),bf3(1,1),nt,nk3,nt,1)
c         ---> form the dot product
         i = 0
         DO 120 i3 = 0,nk3 - 1
            DO 110 i2 = 0,nk2 - 1
               DO 100 i1 = 0,nk1 - 1
                  i = i + 1
                  af3(i,1) = af3(i,1)*af3(i,1)
                  bf3(i,1) = 0.
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
c         ---> back fft
         CALL cfft(af3(1,1),bf3(1,1),nt,nk1,nk1,-1)
         CALL cfft(af3(1,1),bf3(1,1),nt,nk2,nk12,-1)
         CALL cfft(af3(1,1),bf3(1,1),nt,nk3,nt,-1)
         sumis = 0.
         i = 0
         DO 150 i3 = 0,nk3 - 1
            k3 = i3
            IF (k3.GT.mx3) k3 = k3 - nk3
            DO 140 i2 = 0,nk2 - 1
               k2 = i2
               IF (k2.GT.mx2) k2 = k2 - nk2
               DO 130 i1 = 0,nk1 - 1
                  k1 = i1
                  IF (k1.GT.mx1) k1 = k1 - nk1
                  i = i + 1
                  phase = rgphs(k1,k2,k3)
                  id3 = ig(k1,k2,k3)
                  IF (id3.NE.0) THEN
                     sumis = sumis + real(cmplx(af3(i,1),bf3(i,1))*
     +                       conjg(ustep(id3)))*phase*fact
                  END IF
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
         dis(num) = dis(num) + sumis
         IF (film) THEN
c     ----> vacuum part
            DO 220 ivac = 1,nvac
               DO 200 ip = 1,nmzxy
c         ---> create density in the real space
                  i = 0
                  DO 170 i2 = 0,nk2 - 1
                     k2 = i2
                     IF (k2.GT.mx2) k2 = k2 - nk2
                     DO 160 i1 = 0,nk1 - 1
                        k1 = i1
                        IF (k1.GT.mx1) k1 = k1 - nk1
                        i = i + 1
                        id3 = ig(k1,k2,0)
                        IF (id3.EQ.0) THEN
                           af2(i,1) = 0.
                           bf2(i,1) = 0.
                        ELSE
                           id2 = ig2(id3)
                           phase = rgphs(k1,k2,0)
                           IF (id2.EQ.1) THEN
                              af2(i,1) = rhv0(ip,ivac,num,1)
                              bf2(i,1) = 0.
                           ELSE
                              af2(i,1) = real(rhv1(ip,id2-1,ivac,num,1))
                              bf2(i,1) =aimag(rhv1(ip,id2-1,ivac,num,1))
                           END IF
                        END IF
  160                CONTINUE
  170             CONTINUE
                  CALL cfft(af2(1,1),bf2(1,1),nk12,nk1,nk1,1)
                  CALL cfft(af2(1,1),bf2(1,1),nk12,nk2,nk12,1)
c         ---> form dot product
                  i = 0
                  DO 190 i2 = 0,nk2 - 1
                     DO 180 i1 = 0,nk1 - 1
                        i = i + 1
                        af2(i,1) = af2(i,1)*af2(i,1)
                        bf2(i,1) = 0.
  180                CONTINUE
  190             CONTINUE
c         ---> back fft
                  CALL cfft(af2(1,1),bf2(1,1),nk12,nk1,nk1,-1)
                  CALL cfft(af2(1,1),bf2(1,1),nk12,nk2,nk12,-1)
                  disz(npz-ip) = area*af2(1,1)/real(nk12)
  200          CONTINUE
c         ---> beyond warping region
               DO 210 ip = nmzxy + 1,nmz
                  disz(npz-ip) = area*rhv0(ip,ivac,num,1)
  210          CONTINUE
               CALL intgz0(disz,delz,nmz,sumz,tail)
               IF (zrfs .OR. invs) facv = 2.0
               dis(num) = dis(num) + facv*sumz
  220       CONTINUE
         END IF
  230 CONTINUE
      DO 240 num = 1,jspins
         dis(num) = sqrt(dis(num)/vol)*1000.
  240 CONTINUE
      WRITE (6,FMT=8000) iter,dis(1)
      WRITE (16,FMT=8000) iter,dis(1)
 8000 FORMAT (/,'----> distance of  the potential for it=',i3,':',f11.6,
     +       ' mhtr/bohr**3')
      IF (jspins.EQ.2) THEN
         WRITE (6,FMT=8010) iter,dis(2)
         WRITE (16,FMT=8010) iter,dis(2)
 8010    FORMAT (/,'----> distance of spin potential for it=',i3,':',
     +          f11.6,' mhtr/bohr**3')
      END IF

      END SUBROUTINE potdis
      END MODULE m_potdis
