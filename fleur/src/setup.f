      MODULE m_setup
      CONTAINS
      SUBROUTINE setup(
     >     odd,lmaxd,nkptd,jspd,ncvd,natd,ntypd,nlhd,memd,nstd,ntypsd,
     >     nwdd,neigd,nop,nn3d,nn2d,n2d,n3d,k3d,k2d,k1d,kxc3d,kxc2d,
     >     kxc1d,kq3d,kq2d,kq1d,taual,ncst,zelec,phi,nmzxy,nmz,jmtd,
     >     nmzxyd,nmzd,invs2,delz,n_u,l_soc,theta,invs,film,zrfs,l_f,
     >     area,dvac,scale,z1,vol,omtil,l_opti,igrd,jspins,l_ss,ntype,
     >     nwd,nvac,rkmax,gmax,gmaxxc,neq,lmax,amat,bmat,iplot,latnam,
     X     namgrp,volmts,rmt,gw,bbmat,llh,mlh,nmem,clnu,nlh,nsymt,
     <     ntypsy,invarind,invarop,invtab,multab,invsatnr,invsat,
     <     symor,tau,mrot,dx,zatom,jri,namat,
     <     ngopr,nop2,nq2,nq3,ncv,nkpt,lnonsph,nlod,nlo,llod,llo, 
     <     wtkpt,bk,lchg_v,lchange,evac0,el0,ello0,llochg,enmix,
     <     nk3,nk2,nk1,ngz,mx3,mx2,mx1,sk3,sk2,phi2,kv3,kv2,
     <     izmax,izmin,rgphs,ig,igz,ig2,nstr2,nstr,skiplo,
     <     kq1_fft,kq2_fft,kq3_fft,nq2_fft,nq3_fft,kmxq2_fft,
     <     kmxq_fft,kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,
     <     pgft2xy,pgft2yy,pgft2xx,pgft2y,pgft2x,pgfft2,igfft2,
     <     pgfft,igfft,kimax2,kimax,ufft,ustep,zsigma,sigma,sig_b,
     <     d_wgn,ig1,kv1,nstr1,ngopr1,mrot1,tau1,invtab1,multab1,
     <     igfft1,pgfft1,pgft1x,pgft1xx,pgft1xy,pgft1y,pgft1yy)
c
c----------------------------------------
c this routine is called by: fleur.F
c
c setup --+
c         +-- spg2set
c         +-- local_sym -+- ptsym
c         |              +- lhcal -+- gaussp
c         |                        +- gtest
c         |                        +- ylm4
c         +-- strgn1 -+- boxdim
c         |           +- sort
c         |           +- dset
c         |           +- spgrot
c         |           +- unor2or
c         +-- mapatom -+- dotset
c         |            +- dotirl
c         +-- inpeig -- gkptwgt
c         +-- gaunt2 -- grule
c         +-- prp_qfft -+- boxdim
c         |             +- ifft235
c         +-- prp_xcfft -+- boxdim
c         |              +- ifft235
c         +-- stepf
c         +-- convn
c         +-- efield
c----------------------------------------

c
      USE m_localsym
      USE m_rwsymfile
      USE m_spg2set
      USE m_dwigner
      USE m_strgn
      USE m_stepf
      USE m_mapatom
      USE m_writegw
      USE m_convn
      USE m_prpqfft
      USE m_prpxcfft
      USE m_inpeig
      USE m_efield
c-odim
      USE m_od_types, ONLY : od_dim
      USE m_od_mapatom
      USE m_od_chisym
      USE m_od_strgn1
c+odim
      IMPLICIT NONE
C     ..
C     .. Scalars Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,nkptd,jspd,ncvd,natd,ntypd
      INTEGER, INTENT (IN) :: nwdd,nop,nn3d,nn2d,n2d,n3d,k3d,k2d,k1d
      INTEGER, INTENT (IN) :: kxc3d,kxc2d,kxc1d,kq3d,kq2d,kq1d,ntypsd
      INTEGER, INTENT (IN) :: nmzxy,nmz,jmtd,nmzxyd,nmzd,neigd,llod
      INTEGER, INTENT (IN) :: igrd,jspins,nwd,nvac,nlod,nstd,n_u,gw
      INTEGER, INTENT (INOUT) :: ntype,nlhd,memd
      LOGICAL, INTENT (IN) ::film,zrfs,l_f,iplot,invs2,l_opti,l_ss,l_soc
      LOGICAL, INTENT (INOUT) :: invs
      REAL,    INTENT (IN) ::area,dvac,scale,z1,vol,omtil,delz,theta,phi
      REAL,    INTENT (INOUT) :: gmax,gmaxxc
      REAL,    INTENT (OUT) :: zsigma,sigma,sig_b(2)
      CHARACTER*3, INTENT (IN) :: latnam
      CHARACTER*4, INTENT (IN) :: namgrp
      INTEGER, INTENT (OUT) :: nk3,nk2,nk1,ngz,mx3,mx2,mx1
      INTEGER, INTENT (OUT) :: nop2,nq2,nq3,nsymt
      INTEGER, INTENT (OUT) :: kxc1_fft,kxc2_fft,kxc3_fft,kmxxc_fft
      INTEGER, INTENT (OUT) :: nxc3_fft,kimax2,kimax,izmax,izmin
      LOGICAL, INTENT (OUT) :: symor
c-odim
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lmax(ntypd),llo(nlod,ntypd),lnonsph(ntypd)
      INTEGER, INTENT (IN) :: jri(ntypd),nlo(ntypd),ncst(ntypd)
      INTEGER, INTENT (INOUT) :: neq(ntypd)
      REAL,    INTENT (IN) :: zatom(ntypd),dx(ntypd),zelec(nwdd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3)
      REAL,    INTENT (IN) :: taual(3,natd),rmt(ntypd),volmts(ntypd)
      REAL,    INTENT (INOUT) :: rkmax(nwdd)
      REAL,    INTENT (OUT):: tau(3,nop)
      CHARACTER*2, INTENT (IN) :: namat(0:103)
      INTEGER, INTENT (OUT) :: nq2_fft(nwdd),nq3_fft(nwdd)
      INTEGER, INTENT (OUT) :: kmxq2_fft(nwdd),kmxq_fft(nwdd),ncv(ntypd)
      INTEGER, INTENT (OUT) :: kq1_fft(nwdd),kq2_fft(nwdd),kq3_fft(nwdd)
      INTEGER, INTENT (OUT) :: invarind(natd),invarop(natd,nop)
      INTEGER, INTENT (OUT) :: invtab(nop),multab(nop,nop)
      INTEGER, INTENT (OUT) :: invsatnr(natd),invsat(natd),ngopr(natd)
      INTEGER, INTENT (OUT) :: llh(0:nlhd,ntypsd),nkpt(nwdd)
      INTEGER, INTENT (OUT) :: skiplo(ntypd,jspd)
      INTEGER, INTENT (OUT) :: ntypsy(natd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (OUT) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (OUT) :: igfft(0:nn3d-1,2),igfft2(0:nn2d-1,2)
      INTEGER, INTENT (OUT) :: kv2(2,n2d),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (OUT) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (OUT) :: nstr(n3d),nstr2(n2d),ig2(n3d),igz(n3d)
      REAL,    INTENT (OUT) :: sk2(n2d),sk3(n3d),phi2(n2d)
      REAL,    INTENT (OUT) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (OUT) :: pgft2x(0:nn2d-1),pgft2xx(0:nn2d-1)
      REAL,    INTENT (OUT) :: pgft2y(0:nn2d-1),pgft2yy(0:nn2d-1)
      REAL,    INTENT (OUT) :: bbmat(3,3),pgft2xy(0:nn2d-1)
      REAL,    INTENT (OUT) :: pgfft(0:nn3d-1),pgfft2(0:nn2d-1)
      REAL,    INTENT (OUT) :: ufft(0:27*k1d*k2d*k3d-1)
      REAL,    INTENT (OUT) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (OUT) :: bk(3,nkptd,nwdd),wtkpt(nkptd,nwdd)
      REAL,    INTENT (OUT) :: evac0(2,jspd,nwdd),ello0(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: enmix(jspd,nwdd)
      LOGICAL, INTENT (OUT) :: llochg(nlod,ntypd,jspd)
      LOGICAL, INTENT (OUT) :: lchange(0:lmaxd,ntypd,jspd,nwdd)
      LOGICAL, INTENT (OUT) :: lchg_v(2,jspd,nwdd)
      COMPLEX, INTENT (OUT) :: clnu(memd,0:nlhd,ntypsd),ustep(n3d)
      COMPLEX, INTENT (OUT) :: d_wgn(-3:3,-3:3,3,nop)
c-odim
      INTEGER, INTENT (OUT) :: ig1(-odd%k3:odd%k3,-odd%M:odd%M)
      INTEGER, INTENT (OUT) :: kv1(2,odd%n2d),nstr1(odd%n2d)
      INTEGER, INTENT (OUT) :: ngopr1(natd)
      REAL,    INTENT (OUT) :: mrot1(3,3,odd%nop)
      REAL,    INTENT (OUT) :: tau1(3,odd%nop)
      INTEGER, INTENT (OUT) :: invtab1(odd%nop),multab1(odd%nop,odd%nop)
      INTEGER, INTENT (OUT) :: igfft1(0:odd%nn2d-1,2)
      REAL,    INTENT (OUT) :: pgfft1(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgft1x(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgft1xx(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgft1xy(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgft1y(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgft1yy(0:odd%nn2d-1)
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER iwd,nopd,symfh
      REAL rkmaxx
      INTEGER nlhtyp(ntype) 
      CHARACTER(len=1) :: rw
      CHARACTER(len=7) :: symfn
c-odim
      INTEGER ntp1,ii,i,j,n1,n2,na,np1,n
      REAL pps(1:3)
      INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
c+odim
c
      IF (namgrp.eq.'any ') THEN
         rw = 'R'
         symfh = 94 ; symfn = 'sym.out'
         CALL rw_symfile(
     >                   rw,symfh,symfn,nop,bmat,
     X                   mrot,tau,nopd,nop2,symor)
         IF (nop.NE.nopd)  STOP 'setup: nopd =/= nop'
      ELSE
        CALL spg2set(
     >               nop,zrfs,invs,namgrp,latnam,
     <               mrot,tau,nop2,symor)
      ENDIF
      IF (.NOT.odd%d1) THEN
         CALL local_sym(
     >        lmaxd,lmax,nop,mrot,tau,
     >        natd,ntype,neq,amat,bmat,taual,
     X        nlhd,memd,ntypsd,.false.,
     <        nlhtyp,ntypsy,nlh,llh,nmem,mlh,clnu)
         nsymt = ntypsd
c-odim
c        ods%nop = nop
c        ALLOCATE (ods%mrot(1:3,1:3,nop),ods%tau(1:3,nop))
c        ods%mrot(:,:,:) = mrot(:,:,:)
c        ods%tau(:,:) = tau(:,:)
         mrot1(:,:,:) = mrot(:,:,:)
         tau1(:,:) = tau(:,:)
      ELSEIF (odd%d1) THEN
         CALL od_chisym(odd,mrot1,tau1,zrfs,invs,invs2,amat)
         ntp1 = natd
         ALLOCATE (nq1(ntp1),lmx1(ntp1),nlhtp1(ntp1))
         ii = 1
         DO i = 1,ntype
            DO j = 1,neq(i)
               nq1(ii) = 1
               lmx1(ii) = lmax(i)
               ii = ii + 1
            END DO
         END DO
         CALL local_sym(
     >        lmaxd,lmx1,nop,mrot,tau,
     >        natd,ntp1,nq1,amat,bmat,taual,
     <        nlhd,memd,ntypsd,.false.,
     =        nlhtp1,ntypsy,nlh,llh,nmem,mlh,clnu)

         nsymt = ntypsd
         ii = 1
         DO i = 1,ntype
            nlhtyp(i) = nlhtp1(ii)
            ii = ii + neq(i)
         END DO
         DEALLOCATE (lmx1,nlhtp1)
      END IF
c+odim
      IF (n_u.GT.0) THEN
        CALL d_wigner(
     >                nop,mrot,bmat,3,
     <                d_wgn)
      ENDIF
c
c+odim
      IF (.NOT.odd%d1) THEN
         CALL mapatom(
     >             nop,ntypd,natd,ntype,ntypsy,
     >             amat,bmat,mrot,neq,tau,taual,l_f,n_u,
     >             l_soc,theta,phi,
     <             invs,ngopr,invsat,invsatnr,bbmat,
     <             multab,invtab,invarop,invarind)
         ngopr1(1:natd) = ngopr(1:natd)
c        DEALLOCATE ( nq1 )
      ELSE
c-odim
         CALL mapatom(
     >        nop,ntp1,natd,ntp1,ntypsy,
     >        amat,bmat,mrot,nq1,tau,taual,l_f,n_u,
     >        l_soc,theta,phi,
     <        invs,ngopr,invsat,invsatnr,bbmat,
     <        multab,invtab,invarop,invarind)

         CALL od_mapatom(
     >          odd,natd,ntype,taual,neq,mrot1,tau1,amat,bmat,
     <          ngopr1,invtab1,multab1,invsat,invsatnr)
      END IF

c+odim
      IF (film.OR.(namgrp.ne.'any ')) THEN
      CALL strgn1(
     >            k1d,k2d,k3d,n3d,n2d,nop,nn2d,nn3d,natd,
     >            nmzd,jspd,ntypsd,nmzxyd,ntypd,nlhd,jmtd,
     >            nmz,jspins,ntypsy,nmzxy,ntype,nlh,jri,film,zrfs,
     >            invs2,nvac,delz,rmt,dx,namat,zatom,z1,invtab,
     >            igrd,invs,nop2,bmat,gmax,symor,mrot,tau,
     <            mx1,mx2,mx3,nq2,nq3,ngz,nk1,nk2,nk3,neq,
     <            kv2,kv3,sk2,sk3,nstr,nstr2,ig2,igz,ig,
     <            rgphs,izmin,izmax,phi2,
     <            kimax,igfft,pgfft,kimax2,igfft2,pgfft2,
     <            pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy)
c-odim
      IF (odd%d1) THEN
         CALL od_strgn1(
     >     igrd,bmat,zrfs,invs,invs2,odd,
     <     ig1,kv1,nstr1,
     <     igfft1,pgfft1,pgft1x,pgft1xx,pgft1xy,pgft1y,pgft1yy)
      END IF
c+odim
      ELSE
      CALL strgn2(
     >            k1d,k2d,k3d,n3d,n2d,nop,nn2d,nn3d,natd,
     >            nmzd,jspd,ntypsd,nmzxyd,ntypd,nlhd,jmtd,
     >            nmz,jspins,ntypsy,nmzxy,ntype,nlh,jri,film,zrfs,
     >            invs2,nvac,delz,rmt,dx,namat,zatom,z1,invtab,
     >            igrd,invs,bmat,gmax,symor,mrot,tau,
     <            mx1,mx2,mx3,nq2,nq3,ngz,nk1,nk2,nk3,neq,
     <            kv2,kv3,sk2,sk3,nstr,nstr2,ig2,igz,ig,
     <            rgphs,izmin,izmax,
     <            kimax,igfft,pgfft,kimax2,igfft2,pgfft2,
     <            pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy)
      ENDIF
c
      IF (.not.l_opti) THEN
        CALL inpeig(
     >              lmaxd,ntypd,jspd,nkptd,nlod,
     >              ntype,latnam,nwd,jspins,film,nvac,
     >              lmax,amat,bmat,dvac,neq,nlo,odd,
     <              bk,wtkpt,nkpt,ello0,llochg,skiplo,
     <              el0,evac0,lchange,lchg_v,enmix)
      ENDIF
c
c-----> prepare dimensions for charge density fft-box in pwden.f
c
      rkmaxx = 0.0
      DO iwd = 1 , nwd
         CALL  prp_qfft(
     >                  n2d,n3d,kq1d,kq2d,kq3d,
     >                  nq2,nq3,nstr2,nstr,gmax,bmat,sk2,sk3,kv3,l_ss,
     =                  rkmax(iwd),
     <                  kq1_fft(iwd),kq2_fft(iwd),kq3_fft(iwd),
     <                  nq2_fft(iwd),nq3_fft(iwd),
     <                  kmxq2_fft(iwd),kmxq_fft(iwd))
      ENDDO

      IF (gw.ge.1) CALL write_gw(
     >                        ntype,nop,nwd,jspins,natd,
     >                        ncst,neq,lmax,mrot,amat,bmat,rkmax,
     >                        taual,zatom,vol,scale,neigd,lmaxd,
     >                        nlod,llod,nlo,llo)
c
c-----> prepare dimensions for xc fft-box in visxc(g).f
c
      CALL  prp_xcfft(
     >                n3d,kxc1d,kxc2d,kxc3d,
     >                nq3,nstr,gmax,rkmaxx,bmat,sk3,kv3,
     =                gmaxxc,
     <                kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft)
c
      sigma=0.0 ; sig_b(:) = 0.0
      IF (iplot) GOTO 10
c
      CALL stepf(
     >           k1d,k2d,k3d,n3d,natd,ntypd,odd,sk2,n2d,
     >           nq3,ntype,film,neq,kv3,sk3,ig2,omtil,ig,
     >           dvac,area,z1,vol,rmt,taual,bmat,volmts,
     <           ufft,ustep)
c
      CALL convn(
     >           ncvd,ntype,gmax,rmt,
     <           ncv)
c
c--->    set up electric field parameters (if needed) 
c
      CALL efield(
     >            nwdd,ntypd,nstd,ncst,
     >            ntype,nwd,area,neq,zelec,zatom,
     <            zsigma,sigma,sig_b)

  10  CONTINUE

      END SUBROUTINE setup
      END MODULE m_setup
