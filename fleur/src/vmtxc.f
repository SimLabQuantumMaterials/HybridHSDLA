      SUBROUTINE vmtxc(
     >                 jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,nspd,
     >                 clnu,mlh,nmem,llh,nlh,rmsh,ntypsy,jri,natd,neq,
     >                 rho,icorr,total,krla,ntype,jspins,nsymt,
     X                 vr,
     <                 excr)
c     ********************************************************************
c     fit spherical-harmonics expansion of exchange-correlation potential*
c     inside muffint-tin spheres and add it to coulomb potential         *
c                                       c.l.fu and r.podloucky           *
c     ********************************************************************
c     instead of vmtxcor.f: the different exchange-correlation 
c     potentials defined through the key icorr are called through 
c     the driver subroutine vxcall.f, subroutines vectorized
c     ** r.pentcheva 22.01.96
c     *********************************************************
c     angular mesh calculated on speacial gauss-legendre points
c     in order to use orthogonality of lattice harmonics and
c     avoid a least sqare fit
c     ** r.pentcheva 04.03.96
c     *********************************************************
c     In the previous version the parameter statement nspd has
c     been used synonymous to nsp. I have introduced the variable
c     nsp in order to change this. It has cause problems with
c     subroutine: sphglpts.f
c     ** s.bluegel, IFF, 30.July.97
c     *********************************************************

      USE m_lhglpts
      USE m_gaussp
      USE m_xcall, ONLY : vxcall,excall
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,nspd
      INTEGER, INTENT(IN) :: icorr,krla,ntype,jspins,nsymt,natd
      LOGICAL, INTENT(IN) :: total
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (INOUT):: vr(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (OUT)  :: excr(jmtd,0:nlhd,ntypd)
C     ..
C     .. Local Scalars ..
      REAL elh,rr2,r2rho,vlh
      INTEGER jr,js,k,lh,n,nd,i,nsp,nat
C     ..
C     .. Local Arrays ..
      REAL vxc(nspd,jspd),pos0(3),rhoxc(nspd,jspd),
     +     exc(nspd),wt(nspd),rx(3,nspd)
      REAL ylh(nspd,0:nlhd,ntypsd)
C     ..
c     generates nspd points on a sherical shell with radius 1.0
c     angular mesh equidistant in phi, 
c     theta are zeros of the legendre polynomials
c     
      CALL gaussp(
     >            lmaxd,
     <            rx,wt)

      nsp = nspd
c
c     generates the lattice harmonics on the angular mesh 
c     
      CALL lhglpts(
     >             memd,nlhd,ntypsd,lmaxd,nspd,
     >             rx,nsp,
     >             clnu,nmem,mlh,nlh,llh,nsymt,
     <             ylh)
c
c
c     loop over topologically non-equivalent atoms
c
      nat = 1
      DO 200 n = 1,ntype
         nd = ntypsy(nat)
c
c     loop over radial mesh
c
         DO 190 jr = 1,jri(n)
            rr2 = 1.e0/ (rmsh(jr,n)*rmsh(jr,n))
            DO 40 js = 1,jspins
               DO 10 k = 1,nspd
                   rhoxc(k,js) = 0.

   10          CONTINUE

               DO 30 lh = 0,nlh(nd)
                  r2rho = rr2*rho(jr,lh,n,js)
c
c     generate the densities on an angular mesh
c
                  DO 20 k = 1,nsp
                     rhoxc(k,js) = rhoxc(k,js) +
     +                                ylh(k,lh,nd)*r2rho
   20             CONTINUE
   30          CONTINUE
   40       CONTINUE
c
c     calculate the ex.-cor. potential
c
            CALL vxcall
     >                 (6,icorr,krla,jspins,
     >                  nspd,nsp,rhoxc,
     <                  vxc) 
c     ----> now determine the corresponding potential number 
            DO 140 js = 1,jspins
c
c ---> multiplikate vxc with the weights of the k-points    
c
               DO 135 k = 1,nsp
                  vxc(k,js) = vxc(k,js)*wt(k)
 135           CONTINUE   
               DO 130 lh = 0,nlh(nd)
                  vlh = 0
c
c ---> determine the corresponding potential number through gauss integration
c
                  DO 120 k = 1,nsp
                     vlh = vlh + vxc(k,js)*ylh(k,lh,nd)
 120              CONTINUE

c
c ---> add to the given potential

                  vr(jr,lh,n,js) = vr(jr,lh,n,js) + vlh
 130           CONTINUE

  140       CONTINUE
            IF (total) THEN 
c
c     calculate the ex.-cor energy density
c
               CALL excall
     >                      (6,icorr,krla,jspins,
     >                       nspd,nsp,rhoxc,
     <                       exc) 
             END IF
c     ----> now determine the corresponding energy density number 
c
c ---> multiplikate exc with the weights of the k-points    
c
             DO 160 k = 1,nsp
                exc(k) = exc(k)*wt(k)
 160         CONTINUE  
             DO 180 lh = 0,nlh(nd)
                elh = 0
c
c ---> determine the corresponding potential number through gauss integration
c
                DO 170 k = 1,nsp
                   elh = elh + exc(k)*ylh(k,lh,nd)
 170            CONTINUE

c
c ---> add to the given potential

                excr(jr,lh,n) = elh
  180        CONTINUE
   
  190     CONTINUE
          nat = nat + neq(n)
  200  CONTINUE
c

      RETURN
      END
