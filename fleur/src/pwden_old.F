      MODULE m_pwdenold
      PRIVATE 
      PUBLIC pwden_old
      CONTAINS
      SUBROUTINE pwden_old(
     >                     nvd,nbasfcn,nobd,neigd,n3d,jspd,
     >                     k1d,k2d,k3d,
     >                     omtil,mx1,mx2,mx3,nstr,rgphs,ig,
     >                     we,jspin,z,ne,nv,k1,k2,k3,
     <                     qpw)
c     *******************************************************
c     calculates the star function expansion coefficients of
c     the plane wave charge density       c.l.fu
c     *******************************************************

#include"cpp_double.h"
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,nbasfcn,nobd,neigd,n3d,jspd
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,mx1,mx2,mx3
      INTEGER, INTENT (IN) :: jspin,ne,nv
      REAL,    INTENT (IN) :: omtil
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd),nstr(n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: we(nobd)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif
      COMPLEX, INTENT (INOUT):: qpw(n3d,jspd)
C     ..
C     .. Local Scalars ..
      REAL phase,phasep
      INTEGER i,ind,indp,ix1,iy1,iz1,j,n
C     ..
C     .. Local Arrays ..
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL caux1(nobd),caux2(nobd)
#else
      COMPLEX caux1(nobd),caux2(nobd)
#endif
      REAL q1(nvd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC iabs
C     ..
c     ----> g=0 star
c
      DO 10 i = 1,nv
         q1(i) = 0.0
   10 CONTINUE
      DO 30 n = 1,ne
         DO 20 i = 1,nv
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
            q1(i) = q1(i) + we(n)*z(i,n)*z(i,n)
#else
            q1(i)=q1(i)+we(n)*real(z(i,n)*conjg(z(i,n)))
#endif
   20    CONTINUE
   30 CONTINUE
      qpw(1,jspin) = qpw(1,jspin) + ssum(nv,q1,1)/omtil
      DO 60 i = 1,nv
c     ----> g.ne.0 stars
         DO 50 j = 1,i - 1
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
            DO 40 n = 1,ne
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
               caux1(n) = we(n)*z(j,n)*z(i,n)
               caux2(n) = we(n)*z(i,n)*z(j,n)
#else
               caux1(n)=cmplx(we(n),0.0)*z(j,n)*conjg(z(i,n))
               caux2(n)=cmplx(we(n),0.0)*z(i,n)*conjg(z(j,n))
#endif
   40       CONTINUE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
            qpw(ind,jspin) = qpw(ind,jspin) + ssum(ne,caux1,1)*phase
            qpw(indp,jspin) = qpw(indp,jspin) + ssum(ne,caux2,1)*phasep
#else
            qpw(ind,jspin) =qpw(ind,jspin)
     +                  + csum(ne,caux1,1)*cmplx(phase,0.0)
            qpw(indp,jspin)=qpw(indp,jspin)
     +                  + csum(ne,caux2,1)*cmplx(phasep,0.0)
#endif
   50    CONTINUE
   60 CONTINUE

      END SUBROUTINE pwden_old
!--------------------------------------------------------------
      REAL FUNCTION ssum(n,x,incx)

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER incx,n
C     ..
C     .. Array Arguments ..
      REAL x(n)
C     ..
C     .. Local Scalars ..
      INTEGER i
C     ..
      ssum = 0.
      DO i = 1,n,incx
         ssum = ssum + x(i)
      ENDDO
     
      END FUNCTION ssum
!--------------------------------------------------------------
      COMPLEX FUNCTION csum(n,cvec,incx)

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER incx,n
C     ..
C     .. Array Arguments ..
      COMPLEX cvec(n)
C     ..
C     .. Local Scalars ..
      INTEGER i
C     ..
      csum = cmplx(0.0,0.0)
      DO i = 1,n,incx
         csum = csum + cvec(i)
      ENDDO
      END FUNCTION csum
!--------------------------------------------------------------
      END MODULE m_pwdenold
