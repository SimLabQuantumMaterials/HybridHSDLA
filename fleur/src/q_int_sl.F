      MODULE m_qintsl
      CONTAINS
      SUBROUTINE q_int_sl(
     >                   nvd,nbasfcn,neigd,n3d,natd,nop,invtab,
     >                   k1d,k2d,k3d,ntypd,ntype,neq,nq3,
     >                   volmts,volsl,volintsl,taual,symor,tau,mrot,
     >                   omtil,mx1,mx2,mx3,nstr,rgphs,ig,ig2,
     >                   z,ne,nv,k1,k2,k3,rmt,sk3,bmat,kv3,
     >                   nsl,zsl,nmtsl,odi, 
     <                   qintslk)          
c     *******************************************************
c     calculate the charge of the En(k) state 
c     in the interstitial region of each leyer
c                                             Yu.M. Koroteev
c             From pwden_old.F and pwint.F by  c.l.fu
c     *******************************************************
#include"cpp_double.h"
        USE m_constants, ONLY : pimach
        USE m_pwintsl
        USE m_od_types, ONLY : od_inp
c
	IMPLICIT NONE
C
C     .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: nvd,nbasfcn,neigd,n3d,natd,nop
        INTEGER, INTENT (IN) :: k1d,k2d,k3d,ntypd,ntype,mx1,mx2,mx3
        INTEGER, INTENT (IN) :: ne,nv,nq3,nsl
        REAL,    INTENT (IN) :: omtil
        LOGICAL, INTENT (IN) :: symor
C     ..
C     .. Array Arguments ..
        INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd),nstr(n3d)
        INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
        INTEGER, INTENT (IN) :: kv3(3,n3d),neq(ntypd),mrot(3,3,nop)
        INTEGER, INTENT (IN) :: nmtsl(ntypd,natd),invtab(nop)
        REAL,    INTENT (IN) :: volintsl(natd),taual(3,natd),tau(3,nop)
        REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
        REAL,    INTENT (IN) :: zsl(2,natd),volmts(ntypd),volsl(natd) 
        REAL,    INTENT (IN) :: rmt(ntypd),sk3(n3d),bmat(3,3)
        REAL,    INTENT (OUT):: qintslk(nsl,neigd)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Local Scalars ..
      REAL phase,phasep,q1,zsl1,zsl2,qi,volsli,volintsli,tpi
      INTEGER i,ind,indp,ix1,iy1,iz1,j,n,ns
      COMPLEX x
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: stfunint(:,:),z_z(:)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    ALLOCATABLE :: z_h(:,:)
#else
      COMPLEX, ALLOCATABLE :: z_h(:,:)
#endif
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC iabs,conjg,real
      IF (odi%d1) STOP 'well, does not work with 1D. Not
     &                        clear how to define a layer. YM' 
C
      tpi = 2*3.14         !	tpi = 2*pimach(1.0)
c
c     calculate the star function expansion coefficients of
c     the plane wave charge density for each En(k)
c    
c     ----> g=0 star
c
      ALLOCATE ( stfunint(n3d,nsl), z_z(n3d) ) 
c
c  -----> calculate the integrals of star functions over
c                     the layer interstitial
c
      DO i = 1,nsl
         zsl1 = zsl(1,i)
         zsl2 = zsl(2,i)
         volsli = volsl(i)
         volintsli = volintsl(i)
         DO j = 1,nq3
            CALL pwint_sl(
     >                    k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,
     >                    ntype,neq,volmts,taual,zsl1,zsl2,
     >                    volsli,volintsli,symor,tpi,tau,mrot,
     >                    rmt,sk3,bmat,ig2,ig,nmtsl(1,i),
     >                    kv3(1,j),
     <                    x)
            stfunint(j,i) =  x*nstr(j)
       ENDDO  ! over 3D stars
      ENDDO     ! over layers
!
! Here, I reordered the stuff to save memory
!
      DO 30 n = 1,ne
         z_z(:) = cmplx(0.0,0.0)
         q1 = 0.0
         DO 20 i = 1,nv
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
            q1 = q1 + z(i,n)*z(i,n)
#else
            q1 = q1 + real(z(i,n)*conjg(z(i,n)))
#endif
   20    CONTINUE
         z_z(1) = q1/omtil
c
c     ----> g.ne.0 stars
c
      DO 60 i = 1,nv
         DO 50 j = 1,i-1
            ix1 = k1(j) - k1(i)
            iy1 = k2(j) - k2(i)
            iz1 = k3(j) - k3(i)
            IF (iabs(ix1).GT.mx1) GO TO 50
            IF (iabs(iy1).GT.mx2) GO TO 50
            IF (iabs(iz1).GT.mx3) GO TO 50
            ind = ig(ix1,iy1,iz1)
            indp = ig(-ix1,-iy1,-iz1)
            IF (ind.EQ.0 .OR. indp.EQ.0) GO TO 50
            phase = rgphs(ix1,iy1,iz1)/ (nstr(ind)*omtil)
            phasep = rgphs(-ix1,-iy1,-iz1)/ (nstr(indp)*omtil)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
               z_z(ind)  = z_z(ind)  + z(j,n)*z(i,n)*phase
               z_z(indp) = z_z(indp) + z(i,n)*z(j,n)*phasep
#else
               z_z(ind) = z_z(ind) +z(j,n)*conjg(z(i,n))*
     *                         cmplx(phase,0.0)     
               z_z(indp)= z_z(indp)+z(i,n)*conjg(z(j,n))*
     *                         cmplx(phasep,0.0)
#endif
   50    CONTINUE
   60 CONTINUE
c
c ----> calculate a charge in the layer interstitial region of the film
c
      DO i = 1,nsl
         qi = 0.0
         DO j = 1,nq3
            qi = qi + z_z(j)*stfunint(j,i)
         ENDDO            
         qintslk(i,n) = qi 
      ENDDO    ! over layers         
 
  30  ENDDO ! over states

      DEALLOCATE ( stfunint, z_z ) 

      END SUBROUTINE q_int_sl
      END MODULE m_qintsl
