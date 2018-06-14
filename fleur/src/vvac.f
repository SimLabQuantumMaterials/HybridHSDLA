      MODULE m_vvac
c     ****************************************************************
c     calculates the g(2-dim)=0 part of the vacuum coulomb potential *
c     for general symmetry.          c.l.fu, r.podloucky             *
c     ****************************************************************
      CONTAINS
      SUBROUTINE vvac(
     >                nmzd,jspd,n3d,nmz,nvac,dvac,
     >                z1,ng3,zrfs,invs,delz,nstr,bmat,ig2,kv3,
     >                psq,rht,sigma,zsigma,sig_b,area,
     <                vz,rhobar,sig1dh,vz1dh)

      USE m_constants, ONLY : pimach
      USE m_qsf
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nmzd,jspd,n3d
      INTEGER, INTENT (IN) :: nmz,nvac,ng3
      REAL,    INTENT (IN) :: dvac,delz,z1,sigma,zsigma,area
      LOGICAL, INTENT (IN) :: zrfs,invs
      COMPLEX, INTENT (OUT):: rhobar
      REAL,    INTENT (OUT):: sig1dh,vz1dh
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd),bmat(3,3)
      REAL,    INTENT (IN) :: sig_b(2)
      INTEGER, INTENT (IN) :: ig2(n3d),kv3(3,n3d),nstr(n3d)
      COMPLEX, INTENT (IN) :: psq(n3d)
      REAL,    INTENT (OUT):: vz(nmzd,2,jspd)
C     ..
C     .. Local Scalars ..
      COMPLEX sumq,vcons,ci
      REAL bj0,bj1,dh,qzh,fpi,sigmaa(2)
      INTEGER ig3,imz,ivac,ncsh
C     ..
C     .. Local Arrays ..
      REAL f(nmzd),sig(nmzd),vtemp(nmzd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos,sin

      fpi = 4 * pimach()
      ci = cmplx(0.0,1.0)
      
      vz(:,1:nvac,1) = 0.0  ! initialize potential
 
      ! obtain mesh point (ncsh) of charge sheet for external electric field
      !                                                 (X.Nie, IFF, 10/97)
      ncsh = zsigma/delz + 1.01
      sigmaa(1) = (sigma+sig_b(1))/area
      sigmaa(2) = (sigma+sig_b(2))/area
 
      ! g=0 vacuum potential due to neutral charge density
      ! inside slab and zero charge density outside
 
      vcons = -2.*fpi*ci
      dh = z1
      rhobar = -psq(1)
      sumq = (0.0e0,0.0e0)

      DO ig3 = 2,ng3
        IF (ig2(ig3) == 1) THEN           ! select G_|| = 0
          qzh = kv3(3,ig3)*bmat(3,3)*dh
          bj0 = sin(qzh)/qzh
          rhobar = rhobar - psq(ig3)*bj0*nstr(ig3)
          IF (.NOT.(zrfs .OR. invs)) THEN
            bj1 = (sin(qzh) - qzh*cos(qzh)) / (qzh*qzh)
            sumq = sumq + bj1*psq(ig3)*dh*dh
          ENDIF
        ENDIF
      ENDDO

      ivac = 2                        ! lower (ivac=2) vacuum
      IF (nvac.EQ.2) THEN
        vz(1:nmz,ivac,1) = vcons*sumq
      ENDIF

      ! g=0 vacuum potential due to 
      ! negative of rhobar + vacuum (g=0) charge ----> v2(z)
 
      ivac = 1     ! upper vacuum

      CALL qsf(delz,rht(1,ivac,1),sig,nmz,1)

      sig1dh = sig(nmz) - sigmaa(1)  ! need to include contribution from
                                     ! electric field 
      sig(1:nmz) = sig(nmz) - sig(1:nmz)

      CALL qsf(delz,sig,vtemp,nmz,1)

      ! external electric field contribution (X.Nie, IFF, 10/97)
      !                                       corrected 10/99 mw
      DO imz = 1,ncsh
         vz(imz,ivac,1) = -fpi* (vtemp(nmz)-vtemp(imz)) + vz(imz,ivac,1)
     +                    -fpi*(imz-ncsh)*delz*sigmaa(1)
      ENDDO
      DO imz =ncsh+1,nmz
         vz(imz,ivac,1) = -fpi* (vtemp(nmz)-vtemp(imz)) + vz(imz,ivac,1)
      ENDDO

      vz1dh = vz(1,ivac,1)   ! potential on vacuum boundary

      IF (nvac.EQ.1) RETURN

      ivac = 2     ! lower vacuum

      CALL qsf(delz,rht(1,ivac,1),sig,nmz,1)

      f(1:nmz) = sig(1:nmz) - rhobar*dvac + sig1dh 

      CALL qsf(delz,f,vtemp,nmz,1)

      !   external electric field contribution
      !                      WARNING: nvac=2 case NOT verified, beware!
      DO imz = 1,ncsh
         vz(imz,ivac,1) = -fpi* (vtemp(imz)+sig1dh*dvac-
     +                    rhobar*dvac*dvac/2.) + vz1dh + vz(imz,ivac,1)
      ENDDO
      DO imz =ncsh+1,nmz
         vz(imz,ivac,1) = -fpi* (vtemp(imz)+sig1dh*dvac-
     +                    rhobar*dvac*dvac/2.) + vz1dh + vz(imz,ivac,1)
     +                    +fpi*(imz-ncsh)*delz*sigmaa(2)
      ENDDO

      END SUBROUTINE vvac
      END MODULE m_vvac
