      MODULE m_hsohelp
c
c*********************************************************************
c preparation of spin-orbit matrix elements: ahelp, bhelp
c ahelp(i,n,l,m,jspin) =Sum_(G) (conj(c(G,i,jspin)*a(G,n,l,m,jspin))
c bhelp - same a|->b
c Original version replaced by a call to abcof. Maybe not so efficient
c but includes now LO's and could also be used for noco
c                                                        gb`02
c*********************************************************************
c
      CONTAINS
      SUBROUTINE hsohelp(
     >                   jspd,nvd,neigd,ntypd,natd,lmaxd,nop,nbasfcn,
     >                   jspins,nv,nsz,ntype,neq,lmax,ngopr,invtab,
     >                   bmat,omtil,taual,mrot,rmt,rk,k1,k2,k3,bkpt,
     >                   nlod,llod,nlo,llo,l_dulo,nmat,invsat,invsatnr,
     >                   z,us,uds,dus,duds,ulos,dulos,uulon,dulon,ddn,
     >                   lapw_l,bbmat,zso,l_noco,alph,beta,odi,ods,
     >                   nlotot,kveclo,
     <                   ahelp,bhelp,chelp)
c
      USE m_abcof
      USE m_od_types, ONLY : od_inp, od_sym

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,nvd,neigd,ntypd,natd,lmaxd,nop
      INTEGER, INTENT (IN) :: jspins,nv,ntype,nbasfcn,nlod,llod,nmat
      INTEGER, INTENT (IN) :: nlotot
      REAL,    INTENT (IN) :: omtil
      LOGICAL, INTENT (IN) :: l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nsz(jspd),neq(ntypd),lmax(ntypd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),lapw_l(ntypd)
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd),invtab(nop)
      INTEGER, INTENT (IN) :: nlo(ntype),llo(nlod,ntype),kveclo(nlotot)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd)
      REAL,    INTENT (IN) :: rk(nvd),bkpt(3),bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN) :: rmt(ntypd),taual(3,natd),alph(ntypd)
      REAL,    INTENT (IN) ::   us(0:lmaxd,ntypd,jspd),beta(ntypd)
      REAL,    INTENT (IN) ::  dus(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) ::  uds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) :: duds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) ::  ulos(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: dulos(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: uulon(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd,jspd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      COMPLEX, INTENT (INOUT) :: zso(nbasfcn,2*neigd,jspd)
      COMPLEX, INTENT (OUT):: ahelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd)
      COMPLEX, INTENT (OUT):: bhelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd)
      COMPLEX, INTENT (OUT):: chelp(-llod :llod, neigd,nlod,natd, jspd)
#ifdef CPP_INVERSION
      REAL    z(nbasfcn,neigd,jspd)
#else
      COMPLEX z(nbasfcn,neigd,jspd)
#endif
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
!     .. Locals ..
      LOGICAL l_ss
      INTEGER ispin,lmd,l,n,m,na,ie,lm,ll1,nv1(jspd),ngoprl(natd)
      REAL    qss(3)
      INTEGER, ALLOCATABLE :: g1(:,:),g2(:,:),g3(:,:)
      COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:)
!
! turn off the non-collinear part of abcof
!
      l_ss   = .false.
      lmd = lmaxd*(lmaxd+2)
      qss(:) = 0.0
      ngoprl(:) = 1 ! use unrotated coeffs...
!
! some praparations to match array sizes
!
      nv1(1) = nv ; nv1(jspd) = nv
      ALLOCATE ( g1(nvd,jspd),g2(nvd,jspd),g3(nvd,jspd) )
      g1(:,1) = k1(:) ; g1(:,jspd) = k1(:)
      g2(:,1) = k2(:) ; g2(:,jspd) = k2(:)
      g3(:,1) = k3(:) ; g3(:,jspd) = k3(:)

      chelp(:,:,:,:,jspins) = cmplx(0.0,0.0)

      ALLOCATE ( acof(neigd,0:lmd,natd),bcof(neigd,0:lmd,natd) )
      DO ispin = 1, jspins
#if ( defined(CPP_INVERSION) && defined(CPP_SOC) )
         zso(:,1:neigd,ispin) = cmplx(z(:,1:neigd,ispin),0.0)
         CALL abcof(
     >           lmaxd,ntypd,neigd,neigd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngoprl,taual,neq,lmax,rmt,omtil,bmat,
     >           bbmat,bkpt,g1,g2,g3,nv1,nmat,nsz(ispin),zso(1,1,ispin),
     >           us(0,1,ispin),dus(0,1,ispin),uds(0,1,ispin),
     >           duds(0,1,ispin),ddn(0,1,ispin),invsat,invsatnr,
     >           ulos(1,1,ispin),uulon(1,1,ispin),dulon(1,1,ispin),
     >           dulos(1,1,ispin),llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,ispin,alph,beta,qss,kveclo,odi,ods,
     <           acof,bcof,chelp(-llod,1,1,1,ispin))
!
! transfer (a,b)cofs to (a,b)helps used in hsoham
!
         DO ie = 1, neigd
           DO na = 1, natd
             DO l = 1, lmaxd
               ll1 = l*(l+1)
               DO m = -l,l
                lm = ll1 + m
                ahelp(m,l,na,ie,ispin) = (acof(ie,lm,na))
                bhelp(m,l,na,ie,ispin) = (bcof(ie,lm,na))
               ENDDO
             ENDDO
           ENDDO
         ENDDO
         chelp(:,:,:,:,ispin) = (chelp(:,:,:,:,ispin))
#else
         CALL abcof(
     >           lmaxd,ntypd,neigd,neigd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngoprl,taual,neq,lmax,rmt,omtil,bmat,
     >           bbmat,bkpt,g1,g2,g3,nv1,nmat,nsz(ispin),z(1,1,ispin),
     >           us(0,1,ispin),dus(0,1,ispin),uds(0,1,ispin),
     >           duds(0,1,ispin),ddn(0,1,ispin),invsat,invsatnr,
     >           ulos(1,1,ispin),uulon(1,1,ispin),dulon(1,1,ispin),
     >           dulos(1,1,ispin),llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,ispin,alph,beta,qss,kveclo,odi,ods,
     <           acof,bcof,chelp(-llod,1,1,1,ispin))
!
! transfer (a,b)cofs to (a,b)helps used in hsoham
!
         DO ie = 1, neigd
           DO na = 1, natd
             DO l = 1, lmaxd
               ll1 = l*(l+1)
               DO m = -l,l
                lm = ll1 + m
                ahelp(m,l,na,ie,ispin) = (acof(ie,lm,na))
                bhelp(m,l,na,ie,ispin) = (bcof(ie,lm,na))
               ENDDO
             ENDDO
           ENDDO
         ENDDO
#endif
!      write(54,'(6f15.8)')(((chelp(m,ie,1,na,1),m=-1,1),ie=1,5),na=1,2)
!      write(54,'(8f15.8)')(((acof(ie,l,na),l=0,3),ie=1,5),na=1,2)
      ENDDO    ! end of spin loop (ispin)
c
      DEALLOCATE ( acof,bcof,g1,g2,g3 )
      RETURN
      END SUBROUTINE hsohelp
      END MODULE m_hsohelp
