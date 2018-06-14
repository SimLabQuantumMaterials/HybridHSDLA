      MODULE m_intnv
c     ************************************************
c     calculates the integral of charge density 
c     and potential in the unit cell
c     ************************************************
      CONTAINS
      SUBROUTINE int_nv(
     >                  n3d,n2d,nmzxyd,nmzd,jmtd,ntypd,nlhd,ntypsd,
     >                  qpw,vpw_w,
     >                  rhtxy,vxy,
     >                  rht,vz,
     >                  rho,vr,
     >                  area,delz,omtil,rmsh,dx,
     >                  jri,neq,nlh,nstr2,ntypsy,natd,
     >                  ng3,
     >                  nmz,nmzxy,nq2,ntype,nvac,
     >                  zrfs,film,invs,z1,odi,
     <                  result)

#include"cpp_double.h"
      USE m_intgr, ONLY : intgr3,intgz0
      USE m_set,   ONLY : dset
      USE m_od_types, ONLY : od_inp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      REAL area,delz,result,omtil,z1
      INTEGER, INTENT (IN) :: n3d,n2d,nmzxyd,nmzd,natd
      INTEGER, INTENT (IN) :: jmtd,ntypd,nlhd,ntypsd
      INTEGER ng3,
     +        nmz,nmzxy,nq2,ntype,nvac
      LOGICAL film,invs,zrfs
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Array Arguments ..
      COMPLEX qpw(n3d),rhtxy(nmzxyd,odi%n2d-1,2),
     +       vpw_w(n3d),vxy(nmzxyd,odi%n2d-1,2)
      REAL dx(ntypd),rho(jmtd,0:nlhd,ntypd),rht(nmzd,2),
     +     rmsh(jmtd,ntypd),vr(jmtd,0:nlhd,ntypd),vz(nmzd,2)
      INTEGER jri(ntypd),
     +        neq(ntypd),nlh(ntypsd),nstr2(n2d),ntypsy(natd)
C     ..
C     .. Local Scalars ..
      REAL dpdot,facv,tis,tmt,tvac,tvact
      INTEGER i,ip,ivac,j,k2,lh,n,npz,nat
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      REAL dpj(jmtd),dpz(nmzd)
C     ..
C     .. External Functions ..
      COMPLEX CPP_BLAS_cdotc
      EXTERNAL CPP_BLAS_cdotc
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,real
C     ..
C     -----> CALCULATE DENSITY-POTENTIAL INTEGRALS
C
C  ******************* INTERSTITIAL REGION**********************
c
c  -> warping has been moved to vgen and visxc resp. ...gustav
c
      tis = omtil * real( CPP_BLAS_cdotc(ng3,vpw_w,1,qpw,1) )

      WRITE (6,FMT=8020) tis
      WRITE (16,FMT=8020) tis
 8020 FORMAT (/,10x,'interstitial :',t40,f20.10)

      result = result + tis
C
C   ******************M.T. SPHERES*******************
C
      tmt = 0.
      nat = 1
      DO n = 1,ntype
         DO lh = 0,nlh(ntypsy(nat))
            DO j = 1,jri(n)
               dpj(j) = rho(j,lh,n)*vr(j,lh,n)
            ENDDO
            CALL intgr3(dpj,rmsh(1,n),dx(n),jri(n),dpdot)
            tmt = tmt + dpdot*neq(n)
         ENDDO
         nat = nat + neq(n)
      ENDDO 
      WRITE (6,FMT=8030) tmt
      WRITE (16,FMT=8030) tmt
 8030 FORMAT (/,10x,'muffin tin spheres :',t40,f20.10)
      result = result + tmt
C
C *********** VACUUM REGION**************
C
      IF (film .AND. .NOT.odi%d1) THEN
         npz = nmz + 1
         tail = .true.
         IF (zrfs .OR. invs) THEN
            facv = 2.0
         ELSE
            facv = 1.0
         END IF
         tvac = 0.
         tvact = 0.
c     set array dpz to zero
         CALL dset(nmz,0.0,dpz,1)
         DO 150 ivac = 1,nvac
            DO 120 ip = 1,nmz
               dpz(npz-ip) = rht(ip,ivac)*vz(ip,ivac)
C         --->  WARPING REGION
 120        CONTINUE
            DO 140 k2 = 2,nq2
               DO 130 ip = 1,nmzxy
               dpz(npz-ip) = dpz(npz-ip) +
     +                       nstr2(k2)*rhtxy(ip,k2-1,ivac)*
     +                          conjg(vxy(ip,k2-1,ivac))
 130           CONTINUE
 140        CONTINUE
            CALL intgz0(dpz,delz,nmz,tvac,tail)
            tvact = tvact + area*tvac*facv
 150     CONTINUE
         WRITE (6,FMT=8040) tvact
         WRITE (16,FMT=8040) tvact
 8040    FORMAT (/,10x,'vacuum :',t40,f20.10)
         result = result + tvact
       ELSEIF (odi%d1) THEN
c-odim
         npz = nmz +1
         tail = .true.
         tvac = 0.
         tvact = 0.
c     set array dpz to zero
         CALL dset(nmz,0.0,dpz,1)
         DO 125 ip = 1,nmz
            dpz(npz-ip) = (z1+delz*(ip-1))*
     *                    rht(ip,nvac)*vz(ip,nvac)
c          ---> WARPING REGION
 125     CONTINUE
         DO 145 k2 = 2,odi%nq2
            DO 135 ip = 1,nmzxy
               dpz(npz-ip) = dpz(npz-ip)+
     +             (z1+delz*(ip-1))*
     *             rhtxy(ip,k2-1,nvac)*
     *             conjg(vxy(ip,k2-1,nvac))
 135        CONTINUE
 145     CONTINUE
         CALL intgz0(dpz,delz,nmz,tvac,tail)
         tvact = tvact + area*tvac
         WRITE (6,FMT=8041) tvact
         WRITE (16,FMT=8041) tvact
 8041    FORMAT (/,10x,'vacuum :',t40,f20.10)
         result = result + tvact
c+odim
       END IF   

       END SUBROUTINE int_nv
       END MODULE m_intnv
