      MODULE m_alineso
!---------------------------------------------------------------------- 
! Set up SO-hamiltonian for 2nd variation (hsohelp and hsoham) and 
! diagonalize by lapack routines.
! Eigenvalues and vectors (eig_so and zso) are returned 
!----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE alineso(
     > irank,isize,jspd,nvd,neigd,ntypd,natd,lmaxd,nop,nlod,
     > nbasfcn,jspins,ntype,soc_opt,neq,lmax,ngopr,invtab,
     > bmat,omtil,taual,mrot,rmt,
     X rk,bkpt,
     > nlo,llo,llod,
     > rsopp,rsoppd,rsopdp,rsopdpd,nk1,nkpt,l_dulo,
     > rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     > us,dus,uds,duds,ulos,dulos,uulon,dulon,soangl,
     > invsat,invsatnr,ddn,lapw_l,bbmat,l_noco,alph,beta,eonly,
     > nlotot,odi,ods,
     X cp_help,cp_ham,cp_diag,
     < k1,k2,k3,kveclo,nv,wtkpt,ello,nmat,nsize,
     < eig_so,zso)

#include"cpp_double.h"

      USE m_od_types, ONLY : od_inp, od_sym
      USE m_hsohelp
      USE m_hsoham
      USE m_cputime
      USE m_dotir,    ONLY : dotirp
#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : read_eig,read_basis,read_vec,
     +                        read_kptenpar
#endif

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: irank,isize,nlotot
      INTEGER, INTENT (IN) :: jspd,nvd,neigd,ntypd,natd,lmaxd,nop
      INTEGER, INTENT (IN) :: jspins,ntype,nk1,nkpt
      INTEGER, INTENT (IN) :: nlod,nbasfcn,llod
      REAL,    INTENT (IN) :: omtil
      LOGICAL, INTENT (IN) :: l_noco,eonly
      INTEGER, INTENT (OUT):: nmat,nv,nsize
      REAL,    INTENT (OUT):: wtkpt
      REAL,    INTENT (INOUT) :: cp_help,cp_ham,cp_diag
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),lapw_l(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntype),llo(nlod,ntype)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd)
      LOGICAL, INTENT (IN) :: soc_opt(ntypd+2) 
      REAL,    INTENT (INOUT) :: rk(nvd),bkpt(3)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN) :: rmt(ntypd),taual(3,natd),alph(ntypd)
      REAL,    INTENT (IN) ::   us(0:lmaxd,ntypd,jspd),beta(ntypd)
      REAL,    INTENT (IN) ::  dus(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) ::  uds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) :: duds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) ::  ddn(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) :: rsopp  (ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsoppd (ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsopdp (ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsopdpd(ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsoplop (ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsoplopd(ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsopdplo(ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsopplo (ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsoploplop(ntypd,nlod,nlod,2,2)
      REAL,    INTENT (IN) ::  ulos(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: dulos(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: uulon(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd,jspd)
      COMPLEX, INTENT (IN) :: soangl(lmaxd,-lmaxd:lmaxd,2,
     +                               lmaxd,-lmaxd:lmaxd,2)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      COMPLEX, INTENT (OUT) :: zso(nbasfcn,2*neigd,jspd)
      REAL,    INTENT (OUT) :: eig_so(2*neigd),ello(nlod,ntypd,jspd)
      INTEGER, INTENT (OUT) :: k1(nvd),k2(nvd),k3(nvd),kveclo(nlotot)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL      t1,t2,r2
      INTEGER   i,i1,ig,j,jsp,jsp1,k,ne,nn,nn1,nrec,info
      INTEGER   idim_c,idim_r
      CHARACTER vectors 
C     ..
C     .. Local Arrays ..
      INTEGER :: nsz(2)
      REAL    :: bkdu(3),eig(neigd,jspd),s(3)
      REAL    :: epar(0:lmaxd,ntypd),evac(2)
      REAL,   ALLOCATABLE :: rwork(:)
      COMPLEX,ALLOCATABLE :: cwork(:),chelp(:,:,:,:,:)
      COMPLEX,ALLOCATABLE :: ahelp(:,:,:,:,:),bhelp(:,:,:,:,:)
      COMPLEX,ALLOCATABLE :: zhelp(:,:),hso(:,:),hsomtx(:,:,:,:)
#ifdef CPP_INVERSION
      REAL,   ALLOCATABLE :: z(:,:,:)
#else
      COMPLEX,ALLOCATABLE :: z(:,:,:)
#endif
C     ..
C     .. External Subroutines ..
      EXTERNAL CPP_LAPACK_cheev
C     ..
C     .. External Functions ..
      COMPLEX  CPP_BLAS_cdotu,CPP_BLAS_cdotc
      EXTERNAL CPP_BLAS_cdotu,CPP_BLAS_cdotc
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx
C     ..
c     read from eigenvalue and -vector file
c
      ALLOCATE ( z(nbasfcn,neigd,jspd) )

#ifdef CPP_INVERSION 
      z(:,:,:)= 0.  
#else
      z(:,:,:)= CMPLX(0.,0.) 
#endif 
      zso(:,:,:)= CMPLX(0.,0.)

      DO jsp = 1,jspins
        IF (jsp.EQ.1) nrec = nk1 
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
        IF (irank.NE.0) THEN
          IF (jsp.EQ.2) nrec = nk1 + nkpt/isize
        ELSE
          IF (jsp.EQ.2) nrec = nk1 + nkpt
        ENDIF
#else
        IF (jsp.EQ.2) nrec = nk1 + nkpt
#endif
#ifdef CPP_HDF
        CALL read_kptenpar(
     >                     nk1,jsp,
     <                     bkdu,wtkpt,epar,ello(:,:,jsp),
     <                     evac)
        CALL read_eig(
     >                nk1,jsp,
     <                ne,eig(:,jsp))
        CALL read_basis(
     >                  nk1,jsp,nlotot,
     <                  nmat,nv,k1,k2,k3,kveclo)
        CALL read_vec(
     >                nk1,jsp,ne,1,ne,
     <                nmat,z(:,:ne,jsp))
        ! write(*,*) 'process',irank,' reads ',nk1
#else
        READ (66,rec=nrec) epar,evac,ello(:,:,jsp),
     +                     bkdu,wtkpt,ne,nv,nmat,
     +                     eig(:,jsp),k1,k2,k3,kveclo,
     +                     z(:,:,jsp)
        ! write(*,*) 'process',irank,' reads ',nrec
#endif
        bkpt(:) = bkdu(:)
        DO i = 1, nv
          s(1) = bkpt(1) + k1(i)
          s(2) = bkpt(2) + k2(i)
          s(3) = bkpt(3) + k3(i)
          r2 = dotirp(s,s,bbmat)
          rk(i) = sqrt(r2)
        ENDDO

        IF (ne.GT.neigd) THEN
          WRITE (6,'(a13,i4,a8,i4)') 'alineso: ne=',ne,' > neigd=',neigd
          STOP 'alineso: ne > neigd'
        ENDIF
        nsz(jsp) = ne
      ENDDO
c
c set up size of e.v. problem in second variation: nsize
c
      nsize = 0
      DO jsp = 1,jspins
         IF (jspins.EQ.1) THEN
            nsize = 2*nsz(jsp)
            nsz(2) = nsz(1)
         ELSE
            nsize = nsize + nsz(jsp)
         ENDIF
      ENDDO
c
c set up A and B coefficients
c
      ALLOCATE ( ahelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd) )
      ALLOCATE ( bhelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd) )
      ALLOCATE ( chelp(-llod :llod, neigd,nlod,natd ,jspd) )
      CALL cpu_time(t1) 
      CALL hsohelp(
     >             jspd,nvd,neigd,ntypd,natd,lmaxd,nop,nbasfcn,
     >             jspins,nv,nsz,ntype,neq,lmax,ngopr,invtab,
     >             bmat,omtil,taual,mrot,rmt,rk,k1,k2,k3,bkpt,
     >             nlod,llod,nlo,llo,l_dulo,nmat,invsat,invsatnr,
     >             z,us,uds,dus,duds,ulos,dulos,uulon,dulon,ddn,
     >             lapw_l,bbmat,zso,l_noco,alph,beta,odi,ods,
     >             nlotot,kveclo,
     <             ahelp,bhelp,chelp)
      CALL cpu_time(t2)
      cp_help = cp_help + t2 - t1
c
c set up hamilton matrix
c
      ALLOCATE ( hsomtx(2,2,neigd,neigd) )
      CALL hsoham(
     >            jspd,neigd,natd,lmaxd,ntypd,nlod,llod,
     >            ntype,soc_opt,jspins,nsz,lmax,neq,nlo,llo,chelp,
     >            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     >            ahelp,bhelp,rsopp,rsoppd,rsopdp,rsopdpd,soangl,
     <            hsomtx)
      DEALLOCATE ( ahelp,bhelp,chelp )
      CALL cpu_time(t1)
      cp_ham= cp_ham + t1 - t2
c
c add e.v. on diagonal
c
      DO jsp = 1,jspins
         DO i = 1,nsz(jsp)
           hsomtx(jsp,jsp,i,i) = hsomtx(jsp,jsp,i,i) +
     +                           cmplx(eig(i,jsp),0.)
           IF (jspins.EQ.1) THEN
              hsomtx(2,2,i,i) = hsomtx(2,2,i,i) +
     +                           cmplx(eig(i,jsp),0.)
           ENDIF
         ENDDO
      ENDDO
c
c  resort H-matrix 
c
      ALLOCATE ( hso(2*neigd,2*neigd) )
      DO jsp = 1,2
         DO jsp1 = 1,2
            IF (jsp.EQ.1) nn = 0 
            IF (jsp1.EQ.1) nn1 = 0
            IF (jsp.EQ.2) nn = nsz(1)
            IF (jsp1.EQ.2) nn1 = nsz(1)
c
            DO i = 1,nsz(jsp)
              DO j = 1,nsz(jsp1)
                hso(i+nn,j+nn1) = hsomtx(jsp,jsp1,i,j)
              ENDDO
            ENDDO
c
         ENDDO
      ENDDO
      DEALLOCATE ( hsomtx )
c
c diagonalize the hamiltonian using library-routines
c
      idim_c = 4*neigd
      idim_r = 6*neigd

      CALL cpu_time(t1)
      ALLOCATE ( cwork(idim_c),rwork(idim_r) )

      IF (eonly) THEN
        vectors= 'N'
      ELSE
        vectors= 'V'
      ENDIF
      CALL CPP_LAPACK_cheev(vectors,'U',nsize,
     X                      hso,2*neigd,
     <                      eig_so,
     >                      cwork, idim_c, rwork, 
     <                      info)

      IF (info.NE.0) WRITE (6,FMT=8000) info
 8000 FORMAT (' AFTER CPP_LAPACK_cheev: info=',i4)
      CALL cpu_time(t2)
      cp_diag = cp_diag + t2 - t1
      DEALLOCATE ( cwork,rwork )

      IF (.NOT.eonly) THEN
        ALLOCATE ( zhelp(2*neigd,2) )
c
c proj. back to G - space: old eigenvector 'z' to new one 'Z'
c                                 +
c  s      ---    s                | z(G,j,s) ...   z(ig,i,jsp)
c Z (G) = >     z  (G) * C (i,j)  | Z(G,j,s) ... zso(ig,j,jsp)
c  j      ---    i                | C(i,j)   ... hso(i ,j)
c          i                      +
c reorder new e.w.  in 2x2 spin space : zhelp(,1),zhelp(,2)
c
        DO ig = 1,nmat
          DO j = 1,nsize
            DO i1 = 1,2

              jsp = i1
              IF (jspins.EQ.1) jsp = 1
              IF (i1.EQ.1) nn = 0
              IF (i1.EQ.2) nn = nsz(1)
              DO i = 1,nsz(jsp)
                zhelp(i,1) =  hso(i+nn,j)
#ifdef CPP_INVERSION
                zhelp(i,2) =  cmplx(z(ig,i,jsp),0.0)
#else
                zhelp(i,2) =  z(ig,i,jsp)
#endif
              ENDDO
              zso(ig,j,jsp) = CPP_BLAS_cdotc(nsz(jsp),
     +                                   zhelp(1,1),1,zhelp(1,2),1)
            ENDDO  ! jsp
          ENDDO    ! j
        ENDDO      !iG
        DEALLOCATE ( zhelp )
      ENDIF ! (.NOT.eonly)

      DEALLOCATE ( hso,z )
c
      RETURN
      END SUBROUTINE alineso
      END MODULE m_alineso
