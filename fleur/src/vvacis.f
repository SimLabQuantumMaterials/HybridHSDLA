      MODULE m_vvacis
c     **********************************************************
c     g.ne.0 coefficients of vacuum coulomb potential          *
c     due to the interstitial charge density inside slab       *
c                                   c.l.fu, r.podloucky        *
c     **********************************************************
!     modified for thick films to avoid underflows gb`06
!---------------------------------------------------------------
      CONTAINS
      SUBROUTINE vvacis(
     >                  k1d,k2d,k3d,n3d,n2d,jspd,nmzxyd,nvac,nmzxy,ng2,
     >                  zrfs,mx3,z1,delz,bmat,nstr2,nstr,ig,sk2,kv2,
     >                  psq,
     <                  vxy)

      USE m_constants, ONLY : pimach
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,n2d,jspd,nmzxyd
      INTEGER, INTENT (IN) :: nvac,nmzxy,ng2,mx3
      REAL,    INTENT (IN) :: z1,delz
      LOGICAL, INTENT (IN) :: zrfs
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (IN) :: kv2(2,n2d),nstr(n3d),nstr2(n2d)
      REAL,    INTENT (IN) :: bmat(3,3),sk2(n2d)
      COMPLEX, INTENT (IN) :: psq(n3d)
      COMPLEX, INTENT (OUT):: vxy(nmzxyd,n2d-1,2,jspd)
C     ..
C     .. Local Scalars ..
      COMPLEX arg,ci
      REAL dh,g,qz,sign,signz,vcons,z,tpi,minexp,e_m
      REAL arg_r,arg_i
      INTEGER i2d,ig3n,imz,imzxy,ivac,k1,k2,kz,m0,nrec2,nrz,nz
C     ..
C     .. Local Arrays ..
      COMPLEX sumr(2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp
C     ..
      minexp = log(2.0)*MINEXPONENT(tpi)

      tpi = 2 * pimach()
      ci = cmplx(0.0,1.0)

      DO ivac = 1,nvac
         DO i2d = 1,n2d - 1
            DO imzxy = 1,nmzxyd
               vxy(imzxy,i2d,ivac,1) = (0.,0.)
            ENDDO
         ENDDO
      ENDDO
      dh = z1
      m0 = -mx3
      IF (zrfs) m0 = 0
      DO  nrec2 = 2,ng2
         k1 = kv2(1,nrec2)
         k2 = kv2(2,nrec2)
         g = sk2(nrec2)
         IF (-2*dh*g >= minexp ) THEN
           arg_r = exp( - 2*dh*g ) 
         ELSE
           arg_r = 0.0
         ENDIF
         vcons = tpi/g
         DO ivac = 1,nvac
            sumr(ivac) = (0.0,0.0)
            sign = 3. - 2.*ivac
            DO kz = m0,mx3
               ig3n = ig(k1,k2,kz)
c     ----> use only stars within the g_max sphere (oct.97 shz)
               IF (ig3n.ne.0) THEN
                 nz = 1
                 IF (zrfs) nz = nstr(ig3n)/nstr2(nrec2)
                 qz = kz*bmat(3,3)
c     ---> sum over gz-stars
                 DO  nrz = 1,nz
                    signz = 3. - 2.*nrz
                    arg = g + sign*ci*signz*qz
                    arg_i = sign*signz*qz*dh
                    sumr(ivac) = sumr(ivac) + psq(ig3n)*(
     +                   cos(arg_i)*( 1 - arg_r ) + 
     +                ci*sin(arg_i)*( 1 + arg_r ) ) / arg
                 ENDDO  ! nrz 
               ENDIF
            ENDDO       ! kz 
            z = 0 ! moved z1 into above equations gb`06
            DO imz = 1,nmzxy
               IF (-g*z  >= minexp ) THEN
                 e_m = exp( -g*z  )
               ELSE
                 e_m = exp( minexp )
               ENDIF
               vxy(imz,nrec2-1,ivac,1) = vxy(imz,nrec2-1,ivac,1) +
     +                                   vcons*sumr(ivac)*e_m
               z = z + delz
            ENDDO  ! imz 
         ENDDO     ! ivac
      ENDDO        ! nrec2

      END SUBROUTINE vvacis
      END MODULE m_vvacis
