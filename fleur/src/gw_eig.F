      MODULE m_gweig
      CONTAINS
      SUBROUTINE gw_eig(
     >                 lmaxd,ntypd,neigd,nobd,natd,nop,nvd,jspd,
     >                 lmd,nbasfcn,llod,nlod,nlotot,invtab,kveclo,
     >                 ntype,mrot,taual,neq,lmax,lnonsph,rmt,omtil,
     >                 bmat,bbmat,us,dus,uds,duds,ddn,
     >                 invsatnr,invsat,
     >                 ulos,uulon,dulon,dulos,llo,nlo,l_dulo,lapw_l,
     >                 nrec,l_noco,l_ss,jsp,alph,beta,qss,matsize,a2,
     >                 ngopr,jmtd,ntypsy,ntypsd,nlh,nlhd,llh,nmem,memd,
     >                 mlh,clnu,rgphs,ig,nq3,k1d,k2d,k3d,n3d,jri,
     >                 dx,vrtot,vpwtot,vrcou,vpwcou,
     >                 lda_u,vs_mmp,lmaxb,n_u,odi,ods)

#include"cpp_double.h"
      USE m_abcof
      USE m_gwvxc
      USE m_od_types, ONLY : od_inp, od_sym
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,neigd,nobd,natd,nop,nvd,jspd
      INTEGER, INTENT (IN) :: ntype,nbasfcn,llod,nlod,lmd,jri(ntypd)
      INTEGER, INTENT (IN) :: jmtd,ntypsd,nlhd,memd,k1d,k2d,k3d,n3d,nq3
      INTEGER, INTENT (IN) :: n_u,lmaxb ! LDA+U
      REAL,    INTENT (IN) :: omtil
      INTEGER, INTENT (IN) :: jsp,nlotot,matsize,nrec
      LOGICAL, INTENT (IN) :: l_noco,l_ss
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: mrot(3,3,nop),lmax(ntypd),lnonsph(ntypd)
      INTEGER, INTENT (IN) :: lapw_l(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: neq(ntypd),invsatnr(natd),invsat(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),ngopr(natd)
      INTEGER, INTENT (IN) :: ntypsy(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (INOUT) :: kveclo(nlotot)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN) :: taual(3,natd),rmt(ntypd)
      REAL,    INTENT (IN) :: dus(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) ::  us(0:lmaxd,ntypd), uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: ulos(nlod,ntypd),uulon(nlod,ntypd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd),dx(ntypd)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd),qss(3)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: vrtot(jmtd,0:nlhd,ntypd)
      REAL,    INTENT (IN) :: vrcou(jmtd,0:nlhd,ntypd)
      COMPLEX, INTENT (IN) :: vpwtot(n3d),vpwcou(n3d)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
#if ( !defined(CPP_INVERSION) )
      COMPLEX, INTENT (INOUT) :: a2(matsize)
#else
      REAL,    INTENT (IN) :: a2(matsize)
#endif
      INTEGER, INTENT (IN) :: lda_u(ntypd)                          ! LDA+U
      COMPLEX, INTENT (IN) :: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u) !

      INTEGER :: ie,na,n,nn,l,lm,ne,nv(jspd),nmat,i,k,j
      INTEGER :: ngopr_l(natd),invsat_l(natd),nbas
      REAL    :: bkpt(3),ello(nlod,ntypd),wtkpt
      REAL    :: epar(0:lmaxd,ntypd),evac(2)
      INTEGER, ALLOCATABLE :: k1(:,:),k2(:,:),k3(:,:)
      REAL,    ALLOCATABLE :: vxc(:),eig(:)
#if ( !defined(CPP_INVERSION) )
      COMPLEX, ALLOCATABLE :: z(:,:),help(:)
#else
      REAL,    ALLOCATABLE :: z(:,:),help(:)
#endif
      COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)

      character fmt*30,c(0:10)
      data c/'s','p','d','f','g','h','i','j','k','l','m'/
      real r,rr
!
! read in eigenvector of first round
!
      ALLOCATE( z(nbasfcn,neigd),eig(neigd),
     +          k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd) )
      READ (66,rec=nrec) epar,evac,ello(:,:),
     +                   bkpt,wtkpt,ne,nv(jsp),nmat,eig(:),
     +                   k1(:,jsp),k2(:,jsp),k3(:,jsp),kveclo(:),z(:,:)
!
! determine xc-potential matrix elements
!
      ALLOCATE(help(nbasfcn),vxc(ne))
#if ( !defined(CPP_INVERSION) )
      a2 = CONJG(a2)
#endif
      nbas=nv(jsp)+nlotot
      DO ie = 1, ne
        help = 0.0
#ifdef CPP_INVERSION
        CALL CPP_BLAS_sspmv('U',nbas,1.0,a2,z(:,ie),1,
     +                      0.0,help,1)
#else
        CALL CPP_BLAS_chpmv('U',nbas,1.0,a2,z(:,ie),1,
     +                      0.0,help,1)
#endif
        vxc(ie) = eig(ie)-DOT_PRODUCT(z(:,ie),help)
      ENDDO
      WRITE(12,*) vxc(1:ne)
      DEALLOCATE(help,vxc)
!
! calculate abc-coefficients
!
#if ( !defined(CPP_SOC) )
      ALLOCATE( acof(nobd,0:lmd,natd), bcof(nobd,0:lmd,natd),
     +          ccof(-llod:llod,nobd,nlod,natd) )
      CALL abcof(
     >           lmaxd,ntypd,neigd,nobd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngopr,taual,neq,lmax,rmt,omtil,
     >           bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,ne,z,
     >           us,dus,uds,duds,ddn,invsat,invsatnr,
     >           ulos,uulon,dulon,dulos,llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,jsp,alph,beta,qss,kveclo,odi,ods,
     <           acof,bcof,ccof)
      CALL gw_vxc(z,k1(1,jsp),k2(1,jsp),k3(1,jsp),acof,bcof,ccof,! eigenvectors
     >            ne,neigd,nv(jsp),nvd,nbasfcn,                  !      "
     >            neq,ntype,ntypd,natd,jri,jmtd,dx,              ! MT information
     >            lmax,lmaxd,lmd,lnonsph,                        ! l-cutoffs
     >            nlo,nlod,llo,llod,                             ! local orbitals
     >            ntypsy,ntypsd,nlh,nlhd,llh,nmem,memd,mlh,clnu, ! lattice harmonics
     >            rgphs,ig,nq3,k1d,k2d,k3d,n3d,                  ! stars
     >            vrtot,vpwtot,vrcou,vpwcou,                     ! total and coulomb potentials 
     >            lda_u,vs_mmp,lmaxb,n_u,ddn)                    ! LDA+U
      CALL abcrot(ntypd,natd,neigd,lmaxd,lmd,llod,nlod,ntype,neq,
     >           ne,lmax,nlo,llo,nop,ngopr,mrot,invsat,invsatnr,bmat,
     >           odi,
     >           ods,
     X           acof,bcof,ccof)
      na = 0
      DO n = 1, ntype
         DO nn = 1, neq(n)
            na = na + 1
            lm = lmax(n)*(lmax(n)+2) ! lnonsph(n)*(lnonsph(n)+2)
            WRITE(15) lm
            WRITE(15) ((acof(ie,l,na),l=0,lm),ie=1,ne),
     +                ((bcof(ie,l,na),l=0,lm),ie=1,ne),
     &                ((ccof(-llo(i,n):llo(i,n),ie,i,na),i=1,nlo(n)),
     &                  ie=1,ne)
         ENDDO
      ENDDO 
      DEALLOCATE ( acof,bcof,ccof )
      DEALLOCATE ( z,k1,k2,k3,eig )

#else
      STOP 'no SOC + GW'
#endif
      END SUBROUTINE gw_eig
      END MODULE m_gweig
