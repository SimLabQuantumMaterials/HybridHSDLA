      MODULE m_alinemuff
************************************************************************
*                                                                      *
*     eigensystem-solver for moderatly-well converged potentials       *
*     a*z=e*b*z is transformed to h*z'=e*s*z' , whereby                *
*     h=C^T*a*C, s=C^T*b*C and z'=C^(-1)*z, when C is z of the last    *
*     iteration (nv*ne-array)                                          *
*     For ne<<nv the matrixsize is significantly reduced               *
*     aline uses ESSL-calls (use LAPACK's reduc3, tred3, bisect,       *
*     tinvit, trback and rebk3  if no ESSL available):                 *
*     SSPEV:  eigensystem-solver for symmetric, real packes h          *
*             here we have no s-matrix                                 *
*     For all eigenvalues are needed, SSPEV should perform better      *
*     then seclr4 (hope so)                                            *
*                                                     Gustav           *
*                                                                      *
************************************************************************
      CONTAINS
      SUBROUTINE aline_muff(
     >                      nvd,lmaxd,ntypd,neigd,natd,nop,
     >                      lmd,lmplmd,ntype,neq,lmax,rmt,taual,
     >                      bmat,mrot,ngopr,omtil,invsat,
     >                      nb0,jsp,ne,nbasfcn,
     >                      us,dus,uds,duds,
     >                      bkpt,rk,k1,k2,k3,nv,
     X                      z,eig)

#include"cpp_arch.h"
#include"cpp_double.h"

      USE m_hnonmuff
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmd,lmplmd,nvd,lmaxd,ntypd,neigd,natd,nop
      INTEGER, INTENT (IN) :: jsp,nb0,ne,nv,ntype,nbasfcn
      REAL,    INTENT (IN) :: omtil
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: lmax(ntypd),neq(ntypd),ngopr(natd)
      INTEGER, INTENT (IN) :: invsat(natd)
      REAL,    INTENT (IN) :: bkpt(3),rk(nvd),rmt(ntypd),taual(3,natd)
      REAL,    INTENT (IN)  :: bmat(3,3)
      REAL,    INTENT (IN)  :: us(0:lmaxd,ntypd),uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: dus(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (INOUT) :: eig(neigd)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: z(nbasfcn,ne)
#else
      COMPLEX, INTENT (INOUT) :: z(nbasfcn,ne)
#endif
C     ..
C     .. Local Scalars ..
      INTEGER i,info,j,ii
C     ..
C     .. Local Arrays ..
      REAL h(ne*(ne+1)/2),help(3*ne),z1(ne,ne)
C     ..
C     .. External Functions ..
      REAL CPP_BLAS_sdot
      EXTERNAL CPP_BLAS_sdot
C     ..
C     .. External Subroutines ..
      EXTERNAL CPP_LAPACK_ssygv
C     ..
c---> initialize the hamiltonian and overlap matrix
      DO i = 1,ne*(ne+1)/2
            h(i) = 0.0
      END DO

c---> add the diagonal (muffin-tin) terms
      DO i = 1,ne
         ii = (i-1)*i/2 + i
         h(ii) = eig(i)
      END DO

c---> add the off-diagonal (non-muffin-tin) terms
      CALL h_nonmuff(
     >               nvd,lmaxd,ntypd,neigd,natd,nop,
     >               lmd,lmplmd,ntype,neq,lmax,rmt,taual,
     >               bmat,mrot,ngopr,omtil,invsat,
     >               nb0,jsp,z,ne,nbasfcn,
     >               us,dus,uds,duds,
     >               bkpt,rk,k1,k2,k3,nv,
     X               h)

c---> DIAGONALIZE THE HAMILTONIAN USING LIBRARY-ROUTINES

#ifdef CPP_ESSL
c---> ESSL call, IBM AIX
      CALL CPP_LAPACK_sspev (21,
     X                       h,
     <                       eig,z1,
     >                       ne,ne,help,3*ne)
#else
c---> LAPACK call
      CALL CPP_LAPACK_sspev ('V','U',ne,
     X                       h,
     <                       eig,z1,
     >                       ne,help,
     <                       info)
      WRITE (6,FMT=8000) info
 8000 FORMAT (' AFTER CPP_LAPACK_sspev: info=',i4)
#endif


c---> store eigenvectors on array z
      DO i = 1,nv
         DO j = 1,ne
            help(j) = z(i,j)
         END DO
         DO j = 1,ne
            z(i,j) = CPP_BLAS_sdot(ne,help,1,z1(1,j),1)
         END DO
      END DO

      END SUBROUTINE aline_muff
      END MODULE m_alinemuff
