      MODULE m_gwvxc
      CONTAINS
      SUBROUTINE gw_vxc(z,k1,k2,k3,acof,bcof,ccof,                    ! eigenvectors
     >                  nband,neigd,nv,nvd,nbasfcn,                   !      "
     >                  neq,ntype,ntypd,natd,jri,jmtd,dx,             ! MT information
     >                  lmax,lmaxd,lmd,lnonsph,                       ! l-cutoffs
     >                  nlo,nlod,llo,llod,                            ! local orbitals
     >                  ntypsy,ntypsd,nlh,nlhd,llh,nmem,memd,mlh,clnu,! lattice harmonics
     >                  rgphs,ig,nq3,k1d,k2d,k3d,n3d,                 ! stars
     >                  vrtot,vpwtot,vrcou,vpwcou,                    ! total and coulomb potentials 
     >                  lda_u,vs_mmp,lmaxb,n_u,ddn)                   ! LDA+U
C     **********************************************************************
C     * This subroutine calculates expectation values of the vxc-potential *
C     * with the KS wave functions psi(i)                                  *
C     *                                                                    *
C     * vxc(ii) = < psi(i) | vxc | psi(i) >                                *
C     *                                                                    *
C     * and stores them in "vxc" for GW.                                   *
C     *                                                                    *
C     * Christoph Friedrich Oct/2004                                       *
C     **********************************************************************
#include"cpp_double.h"
      USE m_intgr,     ONLY : intgr3
      USE m_constants, ONLY : pimach
      USE m_gaunt,     ONLY : gaunt1
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nvd,nbasfcn,ntypd,natd,jmtd,lmaxd
      INTEGER, INTENT (IN) :: lmd,nlod,llod,ntypsd,nlhd,memd,k1d,k2d
      INTEGER, INTENT (IN) :: k3d,n3d
      INTEGER, INTENT (IN) :: nband,nv,ntype,nq3,n_u,lmaxb
c
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd),jri(ntypd)
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd)
      INTEGER, INTENT (IN) :: lnonsph(ntypd),nlo(ntypd),llo(nlod,ntypd)
      INTEGER, INTENT (IN) :: ntypsy(natd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: dx(ntypd)
      REAL,    INTENT (IN) :: vrtot(jmtd,0:nlhd,ntypd),
     &                        vrcou(jmtd,0:nlhd,ntypd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: acof(neigd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(neigd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,neigd,nlod,natd)
      COMPLEX, INTENT (IN) :: vpwtot(n3d),vpwcou(n3d)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#else
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#endif
      INTEGER, INTENT (IN) :: lda_u(ntypd)
      COMPLEX, INTENT (IN) :: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd)
c
C     .. Local Scalars ..
      INTEGER i_n_u
      INTEGER nlharm,nnbas,nbas,typsym,lm
      INTEGER nbasd,nnbasd
      INTEGER i,l,l1,l2,m,m1,m2,ig1,ig2,j1,j2,ibas1,ibas2
      INTEGER igrid,itype,ilharm,iband,istar,ineq,iatom
      REAL    r
      COMPLEX rc,rr
c
C     .. Local Arrays ..
      INTEGER gg(3)
      REAL    grid(jmtd),rmsh(jmtd),vxc(neigd)
C      COMPLEX :: ustep(n3d)
      COMPLEX vpw(n3d)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX vpwmat(nvd,nvd)
#else
      REAL    vpwmat(nvd,nvd)
#endif
      REAL    vr(jmtd,0:nlhd)
      INTEGER, ALLOCATABLE :: lbas(:)
      REAL,    ALLOCATABLE :: f(:,:,:),integ(:,:,:)
      COMPLEX, ALLOCATABLE :: cof(:),vrmat(:,:)

      nbasd  = 2*(lmaxd+1)+nlod                      ! dimension of MT radial basis functions (l)
      nnbasd = 2*(lmaxd*(lmaxd+2)+1)+nlod*(2*llod+1) ! dimension of MT basis functions (lm)
      allocate(lbas(nbasd),f(jmtd,2,nbasd),cof(nnbasd),
     &  integ(0:nlhd,nbasd,nbasd),vrmat(nnbasd,nnbasd))

C     Rewind unit 14 to first record of MT basis info of current spin.
C     MT basis functions from "radfun" will be ordered as
C     f(::,                   1 : lmax(itype)+1                ) = f
C     f(::,       lmax(itype)+2 : 2*(lmax(itype)+1)            ) = g
C     f(::, 2*(lmax(itype)+1)+1 : 2*(lmax(itype)+1)+nlo(itype) ) = flo
C     (see below)
      m = 2*ntype + sum(lmax(1:ntype)) + sum(nlo(1:ntype))
      DO i=1,m ! rewind to first record of MT basis info of current spin
        BACKSPACE(14)
      ENDDO
      vxc(:)=0

C     Loop over atom types
      iatom=0
      i_n_u=0
      DO itype=1,ntype
        IF(lda_u(itype).ge.0) i_n_u = i_n_u + 1
        nlharm=nlh(ntypsy(sum(neq(:itype-1))+1))
        ! Read in radial basis functions
        READ(14) rmsh(1:jri(itype)) ! radial mesh
        DO l=1,lmax(itype)+1
          READ(14)
     &      f(1:jri(itype),1,l              ),  ! f (large component)
     &      f(1:jri(itype),1,l+lmax(itype)+1),  ! g (large component)
     &      f(1:jri(itype),2,l              ),  ! f (small component)
     &      f(1:jri(itype),2,l+lmax(itype)+1)   ! g (small component)
          lbas(l              )=l-1 ! set pointer to l-quantum number (f)
          lbas(l+lmax(itype)+1)=l-1 ! set pointer to l-quantum number (g)
        ENDDO
        DO i=1,nlo(itype)
          READ(14)
     &      f(1:jri(itype),1,2*(lmax(itype)+1)+i),! flo (large component)
     &      f(1:jri(itype),2,2*(lmax(itype)+1)+i) ! flo (small component)
          lbas(2*(lmax(itype)+1)+i)=llo(i,itype)  ! set pointer to l-quantum number (flo)
        ENDDO
        nbas=2*(lmax(itype)+1)+nlo(itype) ! total number of MT radial basis functions (l)
        ! Calculate vxc = vtot - vcoul
        DO l=0,nlharm
          DO i=1,jri(itype)
            IF(l.eq.0) THEN
              vr(i,0)= -vrcou(i,0,itype)                  !
            ELSE                                          ! vxc = vtot - vcoul
              vr(i,l)=  vrtot(i,l,itype)-vrcou(i,l,itype) !
            ENDIF
          ENDDO
        ENDDO
        ! vr(1:jri(itype),0)=sqrt(2*pimach()) ! for testing
        ! vr(1:jri(itype),1:nlharm)=0         !
          
C       Calculate MT contribution to vxc matrix elements
        ! Precompute auxiliary radial integrals
        DO ilharm=0,nlharm
          DO ibas1=1,nbas
            DO ibas2=1,nbas
              DO igrid=1,jri(itype)
                grid(igrid)=vr(igrid,ilharm)*(
     &            f(igrid,1,ibas1)*f(igrid,1,ibas2)+
     &            f(igrid,2,ibas1)*f(igrid,2,ibas2))
              ENDDO
              CALL intgr3(grid,rmsh,dx(itype),jri(itype),
     &          integ(ilharm,ibas1,ibas2)) ! numerical integration
            ENDDO
          ENDDO
        ENDDO
        
C       Calculate muffin tin contribution to vxc matrix
        vrmat(:,:)=0
        typsym=ntypsy(sum(neq(:itype-1))+1)
        j1=0
        DO ibas1=1,nbas ! loop: left basis function
          l1=lbas(ibas1)
          DO m1=-l1,l1
            j1=j1+1
            j2=0
            DO ibas2=1,nbas ! loop: right basis function
              l2=lbas(ibas2)
              DO 10 m2=-l2,l2
                j2=j2+1
                IF(j2.ne.j1.and.(l1.gt.lnonsph(itype).or.       ! skip off-diagonal elements
     &                           l2.gt.lnonsph(itype))) GOTO 10 ! for l1 or l2 > lnonsph(itype)
c               With the following line instead of the previous one "vxc" will contain the same values as "vxc.old" if there are no local orbitals and lnonsph(*)<lmax(*)
c                if(l1.gt.lnonsph(itype).or.l2.gt.lnonsph(itype)) goto 10
                rr=0
                DO ilharm=0,nlharm ! loop: lattice harmonics of vxc
                  l=llh(ilharm,typsym)
                  DO i=1,nmem(ilharm,typsym)
                    m=mlh(i,ilharm,typsym)
                    rc=clnu(i,ilharm,typsym)*
     &                gaunt1(l1,l,l2,m1,m,m2,lmaxd)
                    rr=rr+integ(ilharm,ibas1,ibas2)*rc
                  ENDDO
                ENDDO
c               LDA+U
                IF(
     &              ibas1.eq.ibas2.and.lda_u(itype).eq.l1.and. ! (LDA+U contrib. only if same radial basis function)
     &              ibas1.le.2*(lmax(itype)+1)                 ! (not for LOs)
     &              ) THEN
                  IF(ibas1.le.lmax(itype)+1) THEN
                    rr=rr+vs_mmp(m2,m1,i_n_u)               ! <f|LDA+U|f> = V
                  ELSE
                    rr=rr+vs_mmp(m2,m1,i_n_u)*ddn(l1,itype) ! <g|LDA+U|g> = V*<g|g>
                  ENDIF
                ENDIF
                rc=cmplx(0,1)**(l2-l1) ! adjusts to a/b/ccof-scaling
                vrmat(j1,j2)=rr*rc
 10           CONTINUE
            ENDDO
          ENDDO
        ENDDO
        nnbas=j1 ! number of MT basis functions (lm)
        
C       Project on eigenvectors
        l=lmax(itype)
        lm=l*(l+2)
        DO ineq=1,neq(itype)
          iatom=iatom+1
          DO iband=1,nband
            cof(   1:lm+1  )=acof(iband,0:lm,iatom)                     !
            cof(lm+2:2*lm+2)=bcof(iband,0:lm,iatom)                     !
            j2=2*lm+2                                                   !
            DO i=1,nlo(itype)                                           ! Store a/b/ccof in the vector cof
              j1=j2+1                                                   !
              j2=j1+2*llo(i,itype)                                      !
              cof(j1:j2)=ccof(-llo(i,itype):llo(i,itype),iband,i,iatom) !
            ENDDO                                                       !
            IF(j2.ne.nnbas) STOP 'gw_vxc: didn''t get NNBAS functions!'
            vxc(iband)=vxc(iband)+sum(conjg(cof(:nnbas))*
     &                 matmul(vrmat(:nnbas,:nnbas),cof(:nnbas)))
          ENDDO
        ENDDO
        
      ENDDO ! End loop over atom types

C     ---------------------------------------------------------------
C     Calculate plane wave contribution
      DO i=1,nq3
        vpw(i)=vpwtot(i)-vpwcou(i) ! calculate vxc=vtot-vcoul
      ENDDO
      ! vpw(:)=ustep(:) !for testing
      
      ! Calculate vxc-matrix,  left basis function (ig1)
      !                       right basis function (ig2)
      i=0
      DO ig1=1,nv
        DO ig2=1,nv
          gg(1)=k1(ig1)-k1(ig2)
          gg(2)=k2(ig1)-k2(ig2)
          gg(3)=k3(ig1)-k3(ig2)
          istar=ig(gg(1),gg(2),gg(3))
          IF(istar.ne.0) THEN
            vpwmat(ig1,ig2)=rgphs(gg(1),gg(2),gg(3))*vpw(istar)
          ELSE
            vpwmat(ig1,ig2)=0
            WRITE(6,'(A,/6I5)') 'Warning: Gi-Gj not in any star:',
     &        k1(ig1),k2(ig1),k3(ig1),k1(ig2),k2(ig2),k3(ig2)
          ENDIF
        ENDDO
      ENDDO

      ! Project on eigenvector
      DO iband=1,nband
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
        vxc(iband)=vxc(iband)+sum(conjg(z(:nv,iband))*
     &             matmul(vpwmat(:nv,:nv),z(:nv,iband)))
#else
        vxc(iband)=vxc(iband)+sum(z(:nv,iband)*
     &             matmul(vpwmat(:nv,:nv),z(:nv,iband)))
#endif
      ENDDO
          
      WRITE(13,*) vxc(1:nband)
      
      deallocate(lbas,f,cof,integ,vrmat)
      
      END SUBROUTINE gw_vxc
      END MODULE m_gwvxc
