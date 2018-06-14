      MODULE m_aline
      CONTAINS
      SUBROUTINE aline(
     >                 nvd,lmaxd,ntypd,neigd,natd,nop,jspd,nwdd,nlod,
     >                 matsize,lmd,lmplmd,ntype,neq,lmax,rmt,taual,
     >                 bmat,mrot,ngopr,omtil,film,invsat,nlotot,
     >                 nrec,jsp,nw,el,lnonsph,us,dus,uds,duds,ddn,
     >                 a,b,k1,k2,k3,rk,nbasfcn,tuu,tud,tdu,tdd,ind,
     >                 invtab,invsatnr,ulos,uulon,dulon,dulos,llod,
     >                 llo,nlo,l_dulo,lapw_l,l_noco,l_ss,alph,beta,
     >                 qss,odi,ods,bbmat,tuulo,tdulo,tuloulo,mlotot,
     >                 mlolotot,ello,cp_sph,cp_nsph,
     <                 cp_seclr4,cp_hssphn,bkpt,nv,z,eig,ne)
*                                                                      *
*     eigensystem-solver for moderatly-well converged potentials       *
*     a*z=e*b*z is transformed to h*z'=e*s*z' , whereby                *
*     h=C^T*a*C, s=C^T*b*C and z'=C^(-1)*z, when C is z of the last    *
*     iteration (nv*ne-array)                                          *
*     For ne<<nv the matrixsize is significantly reduced               *
*     aline uses ESSL-calls (use LAPACK's reduc3, tred3, bisect,       *
*     tinvit, trback and rebk3  if no ESSL available):                 *
*     SSPMV:  matrix-vector multiplication for symmetric matrices      *
*             in packed storage.  	                               *
*     SSYGV:  eigensystem-solver for symmetric, real h and positive    *
*             definite, real, symmetric s using Cholesky-factorisation *
*             tridiagonalisation and a QL-algorithm.                   *
*     For all eigenvalues are needed, DSYGV should perform better      *
*     then seclr4 (hope so)                                            *
*                                                     Gustav           *
*                                                                      *
************************************************************************

#include"cpp_arch.h"
#include"cpp_double.h"

      USE m_abcof
      USE m_hssrwu
      USE m_cputime
      USE m_od_types, ONLY : od_inp, od_sym

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,lmaxd,ntypd,neigd,natd,nop,jspd,nwdd
      INTEGER, INTENT (IN) :: jsp,nw,ntype,lmd,lmplmd,matsize,nrec
      INTEGER, INTENT (IN) :: nlod,nbasfcn,nlotot,llod,mlotot,mlolotot
      REAL,    INTENT (IN) :: omtil
      REAL,    INTENT (INOUT) :: cp_seclr4,cp_hssphn,cp_sph,cp_nsph
      LOGICAL, INTENT (IN) :: film,l_noco,l_ss
      INTEGER, INTENT (OUT):: ne
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: lnonsph(ntypd),lapw_l(ntypd)
      INTEGER, INTENT (IN)  :: k1(nvd),k2(nvd),k3(nvd),mrot(3,3,nop)
      INTEGER, INTENT (IN)  :: lmax(ntypd),neq(ntypd),ngopr(natd)
      INTEGER, INTENT (IN)  :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN)  :: rk(nvd),rmt(ntypd),taual(3,natd)
      REAL,    INTENT (IN)  :: bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN)  :: el(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN)  :: us(0:lmaxd,ntypd),uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: dus(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: ulos(nlod,ntypd),uulon(nlod,ntypd)
      REAL,    INTENT (IN)  :: dulon(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (IN)  :: ddn(0:lmaxd,ntypd),ello(nlod,ntypd)
      REAL,    INTENT (IN)  :: alph(ntypd),beta(ntypd),qss(3)
      REAL,    INTENT (OUT) :: eig(neigd),bkpt(3)
      COMPLEX, INTENT (IN) :: tdd(0:lmplmd,ntypd),tdu(0:lmplmd,ntypd)
      COMPLEX, INTENT (IN) :: tud(0:lmplmd,ntypd),tuu(0:lmplmd,ntypd)
      COMPLEX, INTENT (IN) :: tdulo(0:lmd,-llod:llod,mlotot)
      COMPLEX, INTENT (IN) :: tuulo(0:lmd,-llod:llod,mlotot)
      COMPLEX, INTENT (IN) :: tuloulo(-llod:llod,-llod:llod,mlolotot)
      INTEGER, INTENT (IN) :: ind(0:lmd,0:lmd,ntypd)
      INTEGER, INTENT (IN) :: invtab(nop),invsatnr(natd),invsat(natd)
      INTEGER, INTENT (OUT):: nv(jspd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
#ifdef CPP_INVERSION
      REAL,    INTENT (IN)  :: a(matsize),b(matsize)
      REAL,    INTENT (OUT) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN)  :: a(matsize),b(matsize)
      COMPLEX, INTENT (OUT) :: z(nbasfcn,neigd)
#endif  
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL wtkpt,time1,time2
      INTEGER lhelp
      INTEGER i,info,j,nmat
C     ..
C     .. Local Arrays ..
      INTEGER kveclo(nlotot)
      INTEGER, ALLOCATABLE :: kdu1(:),kdu3(:),kdu2(:)
      REAL,    ALLOCATABLE :: edum(:,:),evac(:),eldum(:,:)
      COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)
#ifdef CPP_INVERSION
      REAL,      PARAMETER :: one=1.0, zro=0.0
      REAL,    ALLOCATABLE :: help(:),h(:,:),s(:,:),z1(:,:)
      REAL     CPP_BLAS_sdot
      EXTERNAL CPP_BLAS_sdot,CPP_BLAS_sspmv
#else
      COMPLEX,   PARAMETER :: one=(1.0,0.0), zro=(0.0,0.0)
      COMPLEX, ALLOCATABLE :: help(:),h(:,:),s(:,:),z1(:,:)
      COMPLEX  CPP_BLAS_cdotc
      EXTERNAL CPP_BLAS_cdotc,CPP_BLAS_chpmv
      REAL,    ALLOCATABLE :: rwork(:)
#endif
C     ..
C     .. External Subroutines ..
      EXTERNAL CPP_LAPACK_ssygv
C     ..
c
      ALLOCATE ( kdu1(nvd),kdu3(nvd),kdu2(nvd),evac(2) )
      ALLOCATE ( edum(0:lmaxd,ntypd),eldum(nlod,ntypd) )
      READ (66,rec=nrec) 
     +         edum,evac,eldum,bkpt,wtkpt,ne,nv(jsp),nmat,
     +         eig,kdu1,kdu2,kdu3,kveclo,z
      DEALLOCATE ( kdu1,kdu3,kdu2,evac,edum,eldum)

#ifndef CPP_INVERSION
!     in outeig z is complex conjugated to make it usable for abcof. Here we 
!                       first have to undo this  complex conjugation for the 
      z = conjg(z)    ! multiplication with a and b matrices.
#endif

c
      ALLOCATE ( h(neigd,neigd),s(neigd,neigd) )
      h = zro ; s=zro
#if ( defined(CPP_ESSL) && defined(CPP_INVERSION) )
      ALLOCATE ( help(nbasfcn),z1(neigd,neigd) )
      lhelp=nbasfcn
#else
      IF (nmat.GT.(neigd+2)*neigd) THEN
        ALLOCATE ( help(nmat) )
      ELSE
        ALLOCATE ( help((neigd+2)*neigd) )
      ENDIF
      lhelp= (neigd+2)*neigd
#endif
c
      DO i = 1,ne
#ifdef CPP_INVERSION
         CALL CPP_BLAS_sspmv('U',nmat,one,a,z(1,i),1,zro,
     <                                            help,1)
#else
         CALL CPP_BLAS_chpmv('U',nmat,one,a,z(1,i),1,zro,
     <                                            help,1)
#endif
         DO j = i,ne
#ifdef CPP_INVERSION
            h(j,i) = CPP_BLAS_sdot(nmat,z(1,j),1,help,1)
#else
            h(j,i) = CPP_BLAS_cdotc(nmat,z(1,j),1,help,1)
#endif
         END DO
      END DO

      DO i = 1,ne
#ifdef CPP_INVERSION
         CALL CPP_BLAS_sspmv('U',nmat,one,b,z(1,i),1,zro,
     <                                            help,1)
#else
         CALL CPP_BLAS_chpmv('U',nmat,one,b,z(1,i),1,zro,
     <                                            help,1)
#endif
         DO j = i,ne
#ifdef CPP_INVERSION
            s(j,i) = CPP_BLAS_sdot(nmat,z(1,j),1,help,1)
#else
            s(j,i) = CPP_BLAS_cdotc(nmat,z(1,j),1,help,1)
#endif
         END DO
      END DO

      ALLOCATE ( acof(neigd,0:lmd,natd),bcof(neigd,0:lmd,natd) )
      ALLOCATE ( ccof(-llod:llod,neigd,nlod,natd) ) 
#ifndef CPP_INVERSION
!     conjugate again for use with abcof; finally use cdotc to revert again
      z = conjg(z)
#endif
#ifdef CPP_SOC
      STOP 'no SOC & reduced diagonalization'
#else
      CALL abcof(
     >           lmaxd,ntypd,neigd,neigd,natd,nop,nvd,1,      ! jspd = 1
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngopr,taual,neq,lmax,rmt,omtil,
     >           bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,ne,z,
     >           us,dus,uds,duds,ddn,invsat,invsatnr,
     >           ulos,uulon,dulon,dulos,llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,1,alph,beta,qss,kveclo,odi,ods,  ! ispin = 1
     <           acof,bcof,ccof)
#endif
!
      CALL cpu_time(time2)
      CALL hssr_wu(
     >             nvd,lmaxd,ntypd,neigd,natd,nop,jspd,nwdd,
     >             lmd,lmplmd,ntype,neq,lmax,nbasfcn,invsat,
     >             nlod,llod,jsp,nw,el,lnonsph,ne,ddn,nmat,
     >             tuu,tud,tdu,tdd,tuulo,tdulo,tuloulo,mlotot,
     >             mlolotot,ind,invtab,acof,bcof,ccof,
     X             cp_sph,cp_nsph,h,s)

      DEALLOCATE ( ccof, bcof, acof )
      CALL cpu_time(time1)
      cp_hssphn = cp_hssphn + time1 - time2
c
#ifdef CPP_INVERSION
#ifdef CPP_ESSL
c---> ESSL call, IBM AIX
      CALL CPP_LAPACK_ssygv (1,h,neigd,s,neigd,eig,z1,neigd,ne,help,
     +                       lhelp)
#else
c---> LAPACK call
      CALL CPP_LAPACK_ssygv(1,'V','L',ne,h,neigd,s,neigd,eig,help,lhelp,
     +                      info)
#endif
#else
      ALLOCATE ( rwork(max(1,3*ne-2)) )
      CALL CPP_LAPACK_chegv(1,'V','L',ne,h,neigd,s,neigd,eig,help,lhelp,
     +                      rwork,info)
      DEALLOCATE ( rwork )
#endif
      IF (info /= 0) THEN
        WRITE (6,FMT=8000) info
        IF (i < 0) THEN
          WRITE(6,'(a7,i3,a22)') 'element',info,' has an illegal value'
        ELSEIF (i > ne) THEN
          WRITE(6,'(a2,i3,a22)') 's:',info-ne,' not positive definite'
        ELSE
          WRITE(6,'(a8,i3,a15)') 'argument',info,' not  converged'
        ENDIF
        STOP
      ENDIF
 8000 FORMAT (' AFTER CPP_LAPACK_ssygv: info=',i4)

      CALL cpu_time(time2)
      cp_seclr4 = cp_seclr4 + time2 - time1
c
      DO i = 1,nmat
         DO j = 1,ne
            help(j) = z(i,j)
         END DO
         DO j = 1,ne
#ifdef CPP_INVERSION
#ifdef CPP_ESSL
c--->       for ESSL call
            z(i,j) = CPP_BLAS_sdot(ne,help,1,z1(1,j),1)
#else
c--->       for LAPACK call
            z(i,j) = CPP_BLAS_sdot(ne,help,1,h(1,j),1)
#endif
#else
            z(i,j) = CPP_BLAS_cdotc(ne,help,1,h(1,j),1)
#endif
         END DO
      END DO
      DEALLOCATE ( help,h,s )
#ifdef CPP_ESSL
      DEALLOCATE (z1)
#endif

      END SUBROUTINE aline
      END MODULE m_aline
