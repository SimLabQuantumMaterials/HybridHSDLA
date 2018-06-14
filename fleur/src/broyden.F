      MODULE m_broyden
c################################################################
c     IMIX = 3 : BROYDEN'S FIRST METHOD
c     IMIX = 5 : BROYDEN'S SECOND METHOD
c     IMIX = 7 : GENERALIZED ANDERSEN METHOD
c     sm   : input charge density of iteration m
c            afterwards update rho(m+1)
c     fm   : output minus input charge density of iteration m
c     sm1  : input charge density of iteration m-1
c     fm1   : output minus inputcharge density of iteration m-1
c################################################################
      CONTAINS
      SUBROUTINE broyden(
     >                   z1,jspd,n3d,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,
     >                   nmaph,mapmt,mapvac2,n2d,k1d,k2d,k3d,natd,
     >                   nmap,jspins,l_noco,imix,maxiter,alpha,fm,ufft,
     >                   kimax,igfft,pgfft,nq3,omtil,odi,
     >                   ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >                   film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,
     =                   mit,sm,lpot)

#include"cpp_double.h"
      USE m_metric
      USE m_od_types, ONLY : od_inp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,n3d,ntypd,nmzd,ntypsd,jmtd,natd
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n2d,mmap,mmaph,nmap
      INTEGER, INTENT (IN) :: imix,jspins,maxiter,nmz,nmzxy,nq2,nq3
      INTEGER, INTENT (IN) :: mapmt,mapvac2,ntype,nvac,kimax,n_u
      REAL,    INTENT (IN) :: alpha,area,delz,omtil,z1
      LOGICAL, INTENT (IN) :: l_noco,film,invs,invs2
      INTEGER, INTENT (INOUT) :: mit
      LOGICAL,OPTIONAL,INTENT(IN) :: lpot
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nstr2(n2d),ntypsy(natd)
      INTEGER, INTENT (IN) :: jri(ntypd),neq(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: dx(ntypd),fm(nmap),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
      REAL,    INTENT (INOUT) :: sm(nmap)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER i,it,k,nit,npos,iread,nmaph
      REAL bm,dfivi,fmvm,one,smnorm,vmnorm,zero,alphan
      LOGICAL :: l_pot
C     ..
C     .. Local Arrays ..
      REAL, ALLOCATABLE :: am(:)
      REAL, ALLOCATABLE :: fm1(:),sm1(:),ui(:),um(:),vi(:),vm(:)
C     ..
C     .. External Functions ..
      REAL CPP_BLAS_sdot
      EXTERNAL CPP_BLAS_sdot
C     ..
C     .. External Subroutines ..
      EXTERNAL CPP_BLAS_saxpy,CPP_BLAS_sscal
C     ..
C     .. Data statements ..
      DATA one/1.0e0/,zero/0.0e0/
C
      dfivi = zero
      
      IF (PRESENT(lpot)) THEN
         l_pot=lpot
      ELSE
         l_pot=.false.
      ENDIF

      ALLOCATE (fm1(mmap),sm1(mmap),ui(mmap),um(mmap),vi(mmap),vm(mmap))
      ALLOCATE ( am(maxiter+1) )
c
      IF (mit.ne.1) THEN
c
c     load input charge density (sm1) and  difference of 
c     in and out charge densities (fm1) from previous iteration (m-1)

         REWIND 57
         READ (57) mit,alphan,(fm1(i),i=1,nmap),(sm1(i),i=1,nmap)
         IF ( abs(alpha-alphan) > 0.0001 ) THEN
           WRITE (6,*) 'mixing parameter has been changed; reset'
           WRITE (6,*) 'broyden algorithm or set alpha to',alphan
           STOP 'broyden: mixing parameter (alpha) changed'
         ENDIF
!
!     loop to generate F_m   - F_(m-1)  ... sm1
!                  and rho_m - rho_(m-1) .. fm1
!
         DO k = 1,nmap
            sm1(k) = sm(k) - sm1(k)
            fm1(k) = fm(k) - fm1(k)
         END DO
      END IF
c
c     save F_m and rho_m for next iteration
c
      REWIND 57
      nit = mit +1
      IF (nit > maxiter+1) nit = 1
      WRITE (57) nit,alpha,fm,sm
c
      IF (mit.eq.1) THEN 
c
c     update for rho for mit=1 is straight mixing
c     sm = sm + alpha*fm
c
         CALL CPP_BLAS_saxpy(nmap,alpha,fm,1,sm,1)
      ELSE
c
c     |vi> = w|vi> 
c     loop to generate um : um = sm1 + alpha*fm1 - \sum <fm1|w|vi> ui
c
         DO k = 1,nmap
            um(k) = alpha*fm1(k) + sm1(k)
         END DO
         iread=min(mit-1,maxiter+1)
         DO it = 2,iread
            READ (59,rec=it-1) (ui(i),i=1,nmap),(vi(i),i=1,nmap),dfivi
            am(it) = CPP_BLAS_sdot(nmap,vi,1,fm1,1)
            CALL CPP_BLAS_saxpy(nmap,-am(it),ui,1,um,1)
            write(6,FMT='(5x,"<vi|w|Fm> for it",i2,5x,f10.6)')it,am(it) 
        END DO
c
c     ****************************************
c             broyden's first method
c     ****************************************
c
         IF (imix.eq.3) THEN
c
c     convolute drho(m) with the metric: |fm1> = w|sm1>
c
            CALL metric(
     >           z1,jspd,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,nmaph,
     >           mapmt,mapvac2,jspins,l_noco,sm1,ufft,kimax,igfft,
     >           pgfft,nq3,omtil,k1d,k2d,k3d,n2d,n3d,natd,
     >           ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >           film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,odi,
     <           fm1,l_pot) 
c
c     calculate the norm of sm1 : <sm1|w|sm1>
            smnorm = CPP_BLAS_sdot(nmap,sm1,1,fm1,1)
c
c     loop to generate vm = alpha*sm1  - \sum <ui|w|sm1> vi
c
            DO k = 1,nmap
               vm(k) = alpha*fm1(k) 
            END DO
            DO it = 2,iread
               READ (59,rec=it-1)(ui(i),i=1,nmap),
     +                           (vi(i),i=1,nmap),dfivi
               bm = CPP_BLAS_sdot(nmap,ui,1,fm1,1)
               CALL CPP_BLAS_saxpy(nmap,-bm,vi,1,vm,1)
            write(6,FMT='(5x,"<ui|w|Fm> for it",i2,5x,f10.6)') it, bm 
            END DO
c
c     complete evaluation of vm
c     vmnorm = <um|w|sm1>-<sm1|w|sm1>
c
            vmnorm = CPP_BLAS_sdot(nmap,fm1,1,um,1) - smnorm
c       *   if (vmnorm.lt.tol_10) stop
            CALL CPP_BLAS_sscal(nmap,one/vmnorm,vm,1)
c     write bm(it)
c
         ELSE IF (imix.eq.5) THEN
c
c     ****************************************
c           broyden's second method
c     ****************************************
c
c --> multiply fm1 with metric matrix and store in vm:  w |fm1>  
c
            CALL metric(
     >           z1,jspd,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,nmaph,
     >           mapmt,mapvac2,jspins,l_noco,fm1,ufft,kimax,igfft,pgfft,
     >           nq3,omtil,k1d,k2d,k3d,n2d,n3d,natd,
     >           ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >           film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,odi,
     <           vm,l_pot)

c     
c     calculate the norm of fm1 and normalize vm it: vm = wfm1 / <fm1|w|fm1>
c
            vmnorm = one/CPP_BLAS_sdot(nmap,fm1,1,vm,1)
            CALL CPP_BLAS_sscal(nmap,vmnorm,vm,1)

         ELSE IF (imix.eq.7) THEN
c
c     ****************************************
c           generalized anderson method
c     ****************************************
c     
c     calculate vm = alpha*wfm1 -\sum <fm1|w|vi> <fi1|w|vi><vi|
c     convolute fm1 with the metrik and store in vm
            CALL metric(
     >           z1,jspd,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,nmaph,
     >           mapmt,mapvac2,jspins,l_noco,fm1,ufft,kimax,igfft,pgfft,
     >           nq3,omtil,k1d,k2d,k3d,n2d,n3d,natd,
     >           ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >           film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,odi,
     <           vm,l_pot) 
c
            DO it = 2,iread
               READ (59,rec=it-1) (ui(i),i=1,nmap),
     +                            (vi(i),i=1,nmap),dfivi
               CALL CPP_BLAS_saxpy(nmap,-am(it)*dfivi,vi,1,vm,1)
            END DO
c
            vmnorm = CPP_BLAS_sdot(nmap,fm1,1,vm,1)

c     *   if (vmnorm.lt.tol_10) stop
            CALL CPP_BLAS_sscal(nmap,one/vmnorm,vm,1)
c
c     save dfivi(mit) for next iteration
c
            dfivi = vmnorm
         END IF
c
c     ****************************************
c
c     write um,vm and dfivi on file broyd.?
c
         npos=mit-1
         IF (mit.GT.maxiter+1) THEN
            npos = MOD(mit-2,maxiter)+1
         ENDIF
         WRITE (59,rec=npos) (um(i),i=1,nmap),(vm(i),i=1,nmap),dfivi
c
c     update rho(m+1)
c
c     calculate <fm|w|vm>
c
         fmvm = CPP_BLAS_sdot(nmap,vm,1,fm,1)
         CALL CPP_BLAS_saxpy(nmap,one-fmvm,um,1,sm,1)
      END IF
      mit = mit + 1
      DEALLOCATE ( fm1,sm1,ui,um,vi,vm,am )

      END SUBROUTINE broyden
      END MODULE m_broyden
