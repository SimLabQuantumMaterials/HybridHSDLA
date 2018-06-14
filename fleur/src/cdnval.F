      MODULE m_cdnval
      CONTAINS
      SUBROUTINE cdnval(
     >   ntapwf,npotmatfile,nkpt,jspin,nw,semic,slice,l_noco,
     >   l_ss,l_mperp,cdinf,dos,ndir,vol,volint,volmts,enmix,sk3,vacdos,
     >   integ,layers,izlay,kk,e1s,e2s,nnne,pallst,layerd,neigd,nkptd,
     >   nv2d,ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd,nvd,nmzd,k1d,
     >   k2d,k3d,kq1d,kq2d,kq3d,irank,isize,nbasfcn,nop,n2d,natd,nwdd,
     >   nmzxyd,nspd,lmd,llpd,
     >   kq1_fft,kq2_fft,kq3_fft,nq3_fft,kmxq_fft,
     >   igq_fft,igfft,pgfft,zelec,w,vr,
     >   vrfull,vz,lepr,e1_dos,e2_dos,efermi, 
     >   sig_dos,lchange,lchg_v,l_f,l_geo,bbmat,sfp,invarop,invarind,
     >   multab,invtab,invsat,invsatnr,ntype,film,zrfs,invs,
     >   symor,nvac,nmz,nmzxy,l_soc,jspins,delz,area,omtil,z1,amat,bmat,
     >   taual,pos,nop2,llh,nmem,mlh,clnu,nsymt,ngopr,nlh,ntypsy,mrot,
     >   tau,nstr,nstr2,ig,ig2,rgphs,mx1,mx2,mx3,kv2,kv3,nq3,nq2,rmt,dx,
     >   rmsh,jri,lmax,neq,nlod,llod,nlo,llo,llochg,
     >   lda_u,n_u,skiplo,nlol,lo1l,l_dulo,ulo_der,lapw_l,
     >   tworkf,nstars,nstm,starcoeff,locx,locy,msh,nstd,
     >   odi,ods,sk2,phi2,ncst,zatom,d_wgn,alph,beta,qss,
     X   nrec,n_mmp,cp_time,ello0,el0,evac0,force,
     X   qpw,rhtxy,rho,rht,qvac,qvlay,cdom,cdomvz,cdomvxy,qa21,
     <   chmom,clmom,qis)
c
c     ***********************************************************
c         this subroutin is a modified version of cdnval.F.
c         it calculates a layer charge distribution and an orbital
c         composition for each state in a film.
c         this information is written in the  'ek_orbcomp'  file
c                                    Yu.Koroteev  01.12.2003
c     ***********************************************************
c     flapw7 valence density generator
c                                         c.l.fu
c     zelec(nwdd) used to calculate ne - 6.12.95 r.pentcheva
c
c     changed subroutine to allow parallel writing of vacdos&dosinp
c     used temporary direct access file 84,tmp_dos to store data used
c     in cdninf
c     call of cdninf changed, sympsi is called from cdnval now
c     look for 'ifdef CPP_MPI' -blocks!               d.wortmann 6.5.99
c
c******** ABBREVIATIONS ************************************************
c     nbands   : number of bands in the energy window
c     noccbd   : number of occupied bands
c     slice    : when set to .true. the charge density of a enery range
c                (slice) or a single state is calculated
c     e1s,e2s  : (used if slice) energy range for the slice. if both
c                are set to 0.0 the charge density of the band nr. nnne
c                is calculated
c     pallst   : (used if slice) if set to .true. bands above the
c                Fermi-Energy are taken into account
c     nnne     : (used if slice) number of the band of which the charge
c                density is calculated if e1s and e2s are 0.0
c     kk       : (used if slice) if set to 0 all k-points are used to
c                calculate the charge density of the slice, otherwise
c                only k-points kk is taken into account
c     nslibd   : number of bands in slice
c     ener     : band energy averaged over all bands and k-points,
c                wheighted with the l-like charge of each atom type
c     sqal     : l-like charge of each atom type. sum over all k-points
c                and bands
c***********************************************************************
c
      USE m_types, ONLY : t_orb, t_orbl, t_orblo, t_mt21, t_lo21
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_cotra, ONLY : cotra3
      USE m_radfun
      USE m_radflo
      USE m_rhomt
      USE m_rhonmt
      USE m_rhomtlo
      USE m_rhonmtlo
      USE m_sympsi
      USE m_enpara, ONLY : w_enpara,mix_enpara
      USE m_eparas      ! energy parameters and partial charges
      USE m_qal21       ! off-diagonal part of partial charges
      USE m_abcof
      USE m_topulay
      USE m_nmat        ! calculate density matrix for LDA + U
      USE m_vacden
      USE m_nstm3
      USE m_pwden
      USE m_forcea8
      USE m_forcea12
      USE m_forcea21
      USE m_checkdop    ! check continuity of density on MT radius R
      USE m_int21       ! integrate (spin) off-diagonal radial functions
      USE m_int21lo     ! -"- for u_lo
      USE m_rhomt21     ! calculate (spin) off-diagonal MT-density coeff's
      USE m_rhonmt21    ! -"-                       non-MT-density coeff's
      USE m_cdnmt       ! calculate the density and orbital moments etc. 
      USE m_orbmom      ! coeffd for orbital moments
      USE m_cputime
      USE m_qmtsl       ! These subroutines divide the film into layers
      USE m_qintsl      ! (slabs) and intergate the DOS in these layers
      USE m_slabdim     ! (mt + interstitial)
      USE m_slabgeom    ! (written by Yu.Koroteev, 2003/2004)
      USE m_orbcomp     ! calculate corbital composition (like p_x,p_y,p_z)
      USE m_Ekwritesl   ! and write to file.
      USE m_abcrot2
      USE m_doswrite
      USE m_cylpts
      USE m_cdnread, ONLY : cdn_read0, cdn_read
#ifdef CPP_MPI
      USE m_mpi_col_den ! collect density data from parallel nodes 
      USE m_mpi_col_dos ! collect DOS data from parallel nodes
#endif
      IMPLICIT NONE
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd
      INTEGER, INTENT (IN) :: nvd,nmzd,k1d,k2d,k3d,kq1d,kq2d,kq3d,msh
      INTEGER, INTENT (IN) :: nop,n2d,natd,nwdd,nmzxyd,nspd,lmd,llpd
      INTEGER, INTENT (IN) :: nstd,layerd,neigd,nkptd,nv2d,npotmatfile
      INTEGER, INTENT (IN) :: jspin,ntapwf,nw,ndir,irank,isize
      INTEGER, INTENT (IN) :: kq1_fft,kq2_fft,kq3_fft,nq3_fft,kmxq_fft
      INTEGER, INTENT (IN) :: nkpt,ntype,nvac,nmz,nmzxy,jspins,nop2
      INTEGER, INTENT (IN) :: nsymt,mx1,mx2,mx3,nq3,nq2,n_u
      INTEGER, INTENT (IN) :: layers,kk,nnne,llod,nlod,nbasfcn,lepr
      LOGICAL, INTENT (IN) :: film,zrfs,invs,symor,l_f,l_soc
      LOGICAL, INTENT (IN) :: l_noco,l_ss,l_mperp
      LOGICAL, INTENT (IN) :: dos,semic,slice,cdinf,vacdos,integ,pallst
      REAL, INTENT    (IN) :: e1s,e2s,e1_dos,e2_dos,efermi,sig_dos
      REAL, INTENT    (IN) :: delz,area,omtil,z1,sfp
      REAL, INTENT    (IN) :: vol,volint
      INTEGER, INTENT (INOUT) :: nrec
c+dw
      REAL,    INTENT (IN) :: tworkf,locx(2),locy(2)
      INTEGER, INTENT (IN) :: nstars,nstm
      LOGICAL, INTENT (IN) :: starcoeff
c-dw
C     ..
C     .. Array Arguments ..
c-odim
      REAL, INTENT    (IN) :: sk2(n2d),phi2(n2d)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT(INOUT) :: qpw(n3d,jspd)
      COMPLEX, INTENT(INOUT) :: rhtxy(nmzxyd,odi%n2d-1,2,jspd)
      COMPLEX, INTENT(INOUT) :: cdom(n3d)
      COMPLEX, INTENT(INOUT) :: cdomvz(nmzd,2)
      COMPLEX, INTENT(INOUT) :: cdomvxy(nmzxyd,odi%n2d-1,2)
      COMPLEX, INTENT(INOUT) :: qa21(ntypd)
      INTEGER, INTENT (IN) :: izlay(layerd,2),ngopr(natd),lda_u(ntype)
      INTEGER, INTENT (IN) :: igq_fft(0:kq1d*kq2d*kq3d-1)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,1)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (IN) :: multab(nop,nop),invtab(nop),ncst(ntypd)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd),lapw_l(ntypd)
      INTEGER, INTENT (IN) :: jri(ntypd),kv2(2,n2d),kv3(3,n3d)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),lmax(ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: nlh(ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),ntypsy(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),skiplo(ntypd)
      INTEGER, INTENT (IN) :: nlol(0:llod,ntypd),lo1l(0:llod,ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL, INTENT    (IN) :: zelec(nwdd),tau(3,nop),zatom(ntypd)
      REAL, INTENT    (IN) :: pos(3,natd),dx(ntypd)
      REAL, INTENT    (IN) :: volmts(ntypd),sk3(n3d)
      REAL, INTENT    (IN) :: rmsh(jmtd,ntypd)
      REAL, INTENT    (IN) :: rmt(ntypd),taual(3,natd)
      REAL, INTENT    (IN) :: amat(3,3),bmat(3,3),bbmat(3,3)
      REAL, INTENT    (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL, INTENT    (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL, INTENT    (IN) :: vr(jmtd,ntypd,jspd),vz(nmzd,2)
      REAL, INTENT    (IN) :: vrfull(jmtd,0:nlhd,ntypd)
      REAL, INTENT    (IN) :: w(neigd,nkptd,jspd)
      REAL, INTENT    (IN) :: alph(ntypd),beta(ntypd),qss(3)
      REAL, INTENT    (IN) :: enmix(jspd,nwdd)
      REAL, INTENT   (OUT) :: chmom(ntypd,jspd),clmom(3,ntypd,jspd)
      REAL, INTENT   (OUT) :: qis(neigd,nkptd,jspd)
      REAL, INTENT (INOUT) :: qvac(neigd,2,nkptd,jspd)
      REAL, INTENT (INOUT) :: qvlay(neigd,layerd,2,nkptd,jspd)
      REAL, INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL, INTENT (INOUT) :: rht(nmzd,2,jspd)
      REAL, INTENT (INOUT) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      REAL, INTENT (INOUT) :: evac0(2,jspd,nwdd),ello0(nlod,ntype,jspd)
      REAL, INTENT (INOUT) :: force(3,ntypd,jspd)
      REAL, INTENT (INOUT) :: cp_time(9)
      COMPLEX, INTENT (IN) :: d_wgn(-3:3,-3:3,3,nop)
      COMPLEX, INTENT(INOUT) :: n_mmp(-3:3,-3:3,n_u)
      LOGICAL, INTENT (IN) :: lchange(0:lmaxd,ntypd,jspd,nwdd)
      LOGICAL, INTENT (IN) :: lchg_v(2,jspd,nwdd),llochg(nlod,ntypd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd),l_dulo(nlod,ntypd)
C     ..
C     .. Local Scalars ..
      REAL wk,wronk,time1,time2,sign,emcd_lo,emcd_up
      INTEGER i,ie,iv,ivac,j,k,l,l1,lh,m,n,ilo,isp,
     +        nbands,noded,nodeu,noccbd,nslibd,na,
     +        ikpt,npd,nat,jsp_start,jsp_end,ispin,nrec1
      INTEGER nmat,loplod,skip_t,skip_tt,reclength,reclength2
      INTEGER n_size,i_rec,n_rank,nlotot,ncored,n_start,n_end,noccbd_l
      COMPLEX czero
      LOGICAL l_p0,l_fmpl,l_mcd,l_evp,l_orbcomprot
C     ..
C     .. Local Arrays ..
      INTEGER n_bands(0:neigd),ncore(ntypd)
      REAL    cartk(3),bkpt(3),xp(3,nspd),e_mcd(ntypd,jspins,nstd)
      REAL    ello(nlod,ntypd,jspd),evac(2,jspd)
      REAL    epar(0:lmaxd,ntypd,jspd),evdu(2,jspd)
      REAL    eig(neigd)
      REAL    vz0(2)
      REAL    uuilon(nlod,ntypd),duilon(nlod,ntypd)
      REAL    ulouilopn(nlod,nlod,ntypd)

      INTEGER, PARAMETER :: n2max_nstm3=13
      REAL,    PARAMETER :: emax_nstm3=2.0/27.21
      INTEGER reclength_vw
c+new_sl
        INTEGER nsld,nsl
c
        INTEGER, ALLOCATABLE :: nmtsl(:,:),nslat(:,:)
        REAL,    ALLOCATABLE :: zsl(:,:),volsl(:)
        REAL,    ALLOCATABLE :: volintsl(:)
        REAL,    ALLOCATABLE :: qintsl(:,:,:,:),qmtsl(:,:,:,:)
        REAL,    ALLOCATABLE :: orbcomp(:,:,:,:,:),qmtp(:,:,:,:)
c-new_sl
c-dw
      INTEGER, ALLOCATABLE :: gvac1d(:),gvac2d(:),nv(:),kveclo(:)
      INTEGER, ALLOCATABLE :: jsym(:),ksym(:),k1(:,:),k2(:,:),k3(:,:)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, ALLOCATABLE :: z(:,:)
#else
      REAL,    ALLOCATABLE :: z(:,:)
#endif
      REAL,    ALLOCATABLE :: aclo(:,:,:),acnmt(:,:,:,:,:)
      REAL,    ALLOCATABLE :: bclo(:,:,:),bcnmt(:,:,:,:,:)
      REAL,    ALLOCATABLE :: cclo(:,:,:,:),ccnmt(:,:,:,:,:),we(:)
      REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:)
      REAL,    ALLOCATABLE :: us(:,:,:),uds(:,:,:),ulos(:,:,:)
      REAL,    ALLOCATABLE :: dus(:,:,:),duds(:,:,:),dulos(:,:,:)
      REAL,    ALLOCATABLE :: ddn(:,:,:),uulon(:,:,:),dulon(:,:,:)
      REAL,    ALLOCATABLE :: uloulopn(:,:,:,:),uloulopn21(:,:,:)
      REAL,    ALLOCATABLE :: uu(:,:,:),dd(:,:,:),du(:,:,:)
      REAL,    ALLOCATABLE :: uunmt(:,:,:,:),ddnmt(:,:,:,:)
      REAL,    ALLOCATABLE :: dunmt(:,:,:,:),udnmt(:,:,:,:),sqlo(:,:,:)
      REAL,    ALLOCATABLE :: qal(:,:,:,:),sqal(:,:,:),ener(:,:,:)
      REAL,    ALLOCATABLE :: svac(:,:),pvac(:,:),mcd(:,:,:)
      REAL,    ALLOCATABLE :: enerlo(:,:,:),qmat(:,:,:,:)
      COMPLEX, ALLOCATABLE ::acof(:,:,:,:),bcof(:,:,:,:),ccof(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: acoflo(:,:,:,:),bcoflo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: cveccof(:,:,:,:,:),f_a12(:,:)
      COMPLEX, ALLOCATABLE :: e1cof(:,:,:),e2cof(:,:,:),f_a21(:,:)
      COMPLEX, ALLOCATABLE :: f_b4(:,:),f_b8(:,:)
      COMPLEX, ALLOCATABLE :: aveccof(:,:,:,:),bveccof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uloulop21(:,:,:)
      COMPLEX, ALLOCATABLE :: uunmt21(:,:,:),ddnmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: dunmt21(:,:,:),udnmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: qstars(:,:,:,:),m_mcd(:,:,:,:)
      TYPE (t_orb),  ALLOCATABLE :: orb(:,:,:,:)
      TYPE (t_orbl), ALLOCATABLE :: orbl(:,:,:,:)
      TYPE (t_orblo),ALLOCATABLE :: orblo(:,:,:,:,:)
      TYPE (t_mt21), ALLOCATABLE :: mt21(:,:)
      TYPE (t_lo21), ALLOCATABLE :: lo21(:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg
C     ..
      czero = cmplx(0.0,0.0)
c---> l_fmpl is meant as a switch to to a plot of the full magnet.
c---> density without the atomic sphere approximation for the magnet.
c---> density. It is not completely implemented (lo's missing).
      l_fmpl = .false.
      IF (l_mperp) THEN
c--->    when the off-diag. part of the desinsity matrix, i.e. m_x and
c--->    m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
c--->    is called only once. therefore, several spin loops have been
c--->    added. if l_mperp = F, these loops run only from jspin - jspin.
         jsp_start = 1
         jsp_end   = 2
         ALLOCATE ( mt21(0:lmaxd,ntypd),lo21(nlod,ntypd) )  ! Deallocation at end of subroutine
         ALLOCATE ( uloulopn21(nlod,nlod,ntypd) )
         ALLOCATE ( uloulop21(nlod,nlod,ntypd) )
         ALLOCATE ( qmat(0:3,ntypd,neigd,4) )
         IF (l_fmpl) THEN
            ALLOCATE ( uunmt21((lmaxd+1)**2,nlhd,ntypd) )
            ALLOCATE ( ddnmt21((lmaxd+1)**2,nlhd,ntypd) )
            ALLOCATE ( dunmt21((lmaxd+1)**2,nlhd,ntypd) )
            ALLOCATE ( udnmt21((lmaxd+1)**2,nlhd,ntypd) )
         ELSE
            ALLOCATE ( uunmt21(1,1,1),ddnmt21(1,1,1) )
            ALLOCATE ( dunmt21(1,1,1),udnmt21(1,1,1) )
         ENDIF
      ELSE
         jsp_start = jspin
         jsp_end   = jspin
         ALLOCATE ( mt21(1,1),lo21(1,1),uunmt21(1,1,1) )
         ALLOCATE ( ddnmt21(1,1,1),dunmt21(1,1,1),udnmt21(1,1,1) )
         ALLOCATE ( uloulopn21(1,1,1),uloulop21(1,1,1),qmat(1,1,1,1) )
      ENDIF
c
c---> if l_mperp = F, these variables are only needed for one spin 
c---> at a time, otherwise for both spins:
c
      ALLOCATE ( f(jmtd,2,0:lmaxd,jsp_start:jsp_end) )      ! Deallocation before mpi_col_den
      ALLOCATE ( g(jmtd,2,0:lmaxd,jsp_start:jsp_end) )
      ALLOCATE (   us(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE (  uds(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE (  dus(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( duds(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( ddn(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd) )   ! Deallocation at end of subroutine
      ALLOCATE ( jsym(neigd),ksym(neigd) )
      ALLOCATE ( gvac1d(nv2d),gvac2d(nv2d),nv(jspd) )
      ALLOCATE (  ulos(nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( dulos(nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( uulon(nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( dulon(nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( uloulopn(nlod,nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( uu(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( dd(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( du(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( uunmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( ddnmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( dunmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( udnmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( qal(0:3,ntypd,neigd,jsp_start:jsp_end) )
      ALLOCATE ( sqal(0:3,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( ener(0:3,ntypd,jsp_start:jsp_end) )
      ALLOCATE (   sqlo(nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( enerlo(nlod,ntypd,jsp_start:jsp_end) )
      ALLOCATE ( svac(2,jsp_start:jsp_end) )
      ALLOCATE ( pvac(2,jsp_start:jsp_end) )
      ALLOCATE ( qstars(nstars,neigd,layerd,2) )
!
! --> Initializations
!
      uu(:,:,:) = 0.0 ; dd(:,:,:) = 0.0 ; du(:,:,:) = 0.0
      IF (l_mperp) THEN
         mt21(:,:)%uu = czero ; mt21(:,:)%ud = czero 
         mt21(:,:)%du = czero ; mt21(:,:)%dd = czero 
         lo21(:,:)%uulo = czero ; lo21(:,:)%ulou = czero
         lo21(:,:)%dulo = czero ; lo21(:,:)%ulod = czero
         uloulop21(:,:,:) = czero
      ENDIF
      uunmt(:,:,:,:) = 0.0 ; ddnmt(:,:,:,:) = 0.0
      udnmt(:,:,:,:) = 0.0 ; dunmt(:,:,:,:) = 0.0
      IF (l_fmpl) THEN
         IF (.not.l_mperp) STOP 'for fmpl set l_mperp=T!'
         uunmt21(:,:,:) = czero ; ddnmt21(:,:,:) = czero
         udnmt21(:,:,:) = czero ; dunmt21(:,:,:) = czero
      ENDIF
      svac(:,:) = 0.0 ; pvac(:,:) = 0.0
      sqal(:,:,:) = 0.0 ; ener(:,:,:) = 0.0
c+soc
      IF (l_soc) THEN
        ALLOCATE ( orb(0:lmaxd,-lmaxd:lmaxd,ntypd,jsp_start:jsp_end) )   ! Deallocation at end of subroutine
        ALLOCATE ( orbl(nlod,-llod:llod,ntypd,jsp_start:jsp_end)     )
        ALLOCATE ( orblo(nlod,nlod,-llod:llod,ntypd,jsp_start:jsp_end))
        orb(:,:,:,:)%uu = 0.0 ; orb(:,:,:,:)%dd = 0.0
        orb(:,:,:,:)%uum = czero ; orb(:,:,:,:)%uup = czero
        orb(:,:,:,:)%ddm = czero ; orb(:,:,:,:)%ddp = czero
        orbl(:,:,:,:)%uulo = 0.0 ; orbl(:,:,:,:)%dulo = 0.0
        orbl(:,:,:,:)%uulom = czero ; orbl(:,:,:,:)%uulop = czero
        orbl(:,:,:,:)%dulom = czero ; orbl(:,:,:,:)%dulop = czero
        orblo(:,:,:,:,:)%z = 0.0
        orblo(:,:,:,:,:)%p = czero ; orblo(:,:,:,:,:)%m = czero
      ELSE
        ALLOCATE ( orb(0:0,-lmaxd:-lmaxd,1,jsp_start:jsp_end) )
        ALLOCATE ( orbl(1,-llod:-llod,1,jsp_start:jsp_end) )
        ALLOCATE ( orblo(1,1,-llod:-llod,1,jsp_start:jsp_end) )
      ENDIF
c+for
      IF (l_f) THEN
         ALLOCATE ( f_a12(3,ntypd),f_a21(3,ntypd) )           ! Deallocation at end of subroutine
         ALLOCATE ( f_b4(3,ntypd),f_b8(3,ntypd) )
         f_b4(:,:) = czero  ; f_a12(:,:) = czero
         f_b8(:,:) = czero  ; f_a21(:,:) = czero
      ELSE
         ALLOCATE ( f_b8(1,1) )
      ENDIF
c
      l_mcd = .false.
      INQUIRE (file='mcd_inp',exist=l_mcd)
      IF (l_mcd) THEN
        OPEN (23,file='mcd_inp',STATUS='old',FORM='formatted')
        READ (23,*) emcd_lo,emcd_up
        CLOSE (23)
        ALLOCATE ( m_mcd(nstd,(3+1)**2,3*ntypd,2) )           ! Deallocation at end of subroutine
        ALLOCATE ( mcd(3*ntypd,nstd,neigd) )
        IF (.not.dos) WRITE (*,*) 'For mcd-spectra set dos=T!'
      ELSE
        ALLOCATE ( m_mcd(1,1,1,1),mcd(1,1,1) )
      ENDIF
c
c     determine the total number of lo's : nlotot
c
      nlotot = 0
      DO n = 1, ntype
        DO l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        ENDDO
      ENDDO
      ALLOCATE ( kveclo(nlotot) )

      l_p0 = .false.
      IF (irank.EQ.0) l_p0 = .true.
      IF (l_p0) WRITE (6,FMT=8000) jspin
      IF (l_p0) WRITE (16,FMT=8000) jspin
 8000 FORMAT (/,/,10x,'valence density: spin=',i2)

      CALL cdn_read0(
     >               lmaxd,ntypd,nlod,neigd,jspd,
     >               irank,isize,jspin,jsp_start,jsp_end,
     >               l_noco,nrec,ntapwf,
     <               ello,evac,epar,bkpt,wk,n_bands,nrec1,n_size)

c+lo
c---> if local orbitals are used, the eigenvector has a higher 
c---> dimension then nvd
      loplod  = nlod*(nlod+1)/2
      ALLOCATE ( aclo(nlod,ntypd,jsp_start:jsp_end),                     ! Deallocated at end of subroutine
     +           bclo(nlod,ntypd,jsp_start:jsp_end),
     +           cclo(nlod,nlod,ntypd,jsp_start:jsp_end),
     +           acnmt(0:lmaxd,nlod,nlhd,ntypd,jsp_start:jsp_end), 
     +           bcnmt(0:lmaxd,nlod,nlhd,ntypd,jsp_start:jsp_end), 
     +           ccnmt(nlod,nlod,nlhd,ntypd,jsp_start:jsp_end) )
      aclo(:,:,:) = 0.0 ; bclo(:,:,:) = 0.0 ; ccnmt(:,:,:,:,:) = 0.0
      acnmt(:,:,:,:,:)=0.0 ; bcnmt(:,:,:,:,:)=0.0 ; cclo(:,:,:,:)=0.0
      skip_tt = 0
      DO n = 1,ntype
        skip_tt = skip_tt + skiplo(n)*neq(n)
      ENDDO
      IF (l_soc.OR.l_noco)  skip_tt = 2 * skip_tt
c-lo
c---> set up l-dependent m.t. wavefunctions
      na = 1
      ncored = 0
      ALLOCATE ( flo(jmtd,2,nlod,jspd) )
      DO 40 n = 1,ntype
         IF (cdinf.AND.l_p0) WRITE (6,FMT=8001) n
         DO 30 l = 0,lmax(n)
            DO ispin = jsp_start,jsp_end
               CALL radfun(
     >              l,epar(l,n,ispin),vr(1,n,ispin),jri(n),rmsh(1,n),
     >              dx(n),jmtd,
     <              f(1,1,l,ispin),g(1,1,l,ispin),us(l,n,ispin),
     <              dus(l,n,ispin),uds(l,n,ispin),duds(l,n,ispin),
     <              ddn(l,n,ispin),nodeu,noded,wronk)
               IF (cdinf.AND.l_p0) WRITE (6,FMT=8002) l,
     +              epar(l,n,ispin),us(l,n,ispin),dus(l,n,ispin),nodeu,
     +              uds(l,n,ispin),duds(l,n,ispin),noded,ddn(l,n,ispin),
     +              wronk
            ENDDO
            IF (l_mperp) THEN
               CALL int_21(
     >                     f,g,rmsh(1,n),dx(n),jri(n),jmtd,lmaxd,jspd,l,
     <                     mt21(l,n)%uun,mt21(l,n)%udn,
     <                     mt21(l,n)%dun,mt21(l,n)%ddn)
            ENDIF
   30    CONTINUE
         IF (l_mcd) THEN
            CALL mcd_init(
     >                    lmaxd,jspins,jmtd,msh,ntypd,nstd,
     >                    jri,lmax,dx,zatom,rmsh,ncst,
     >                    vr,g,f,emcd_up,emcd_lo,n,jspin,
     <                    ncore,e_mcd,m_mcd)
            ncored = max(ncore(n),ncored)
         ENDIF
c
c--->   generate the extra wavefunctions for the local orbitals,
c--->   if there are any.
c
         DO ilo = 1, nlo(n)
            DO ispin = jsp_start,jsp_end
              CALL radflo(
     >                    ntypd,nlod,jspd,jmtd,lmaxd,n,ispin,
     >                    ello(1,1,ispin),vr(1,n,ispin),
     >                    jri(n),rmsh(1,n),dx(n),f(1,1,0,ispin),
     >                    g(1,1,0,ispin),llo,nlo,l_dulo(1,n),irank,
     >                    ulo_der,  
     <                    ulos(1,1,ispin),dulos(1,1,ispin),
     <                    uulon(1,1,ispin),dulon(1,1,ispin),
     <                    uloulopn(1,1,1,ispin),
     <                    uuilon,duilon,ulouilopn,  
     <                    flo(1,1,1,ispin))
            ENDDO
            IF (l_mperp) THEN
               CALL int_21lo(
     >                       f,g,rmsh(1,n),dx(n),jri(n),jmtd,lmaxd,jspd,
     >                       flo,llo(1,n),ilo,nlo(n),nlod,
     <                       lo21(ilo,n)%uulon,lo21(ilo,n)%dulon,
     <                       lo21(ilo,n)%uloun,lo21(ilo,n)%ulodn,
     <                       uloulopn21(1,1,n))
            ENDIF
         ENDDO

         na = na + neq(n)
   40 CONTINUE
      DEALLOCATE (flo)
 8001 FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',/,
     +       t32,'radial function',t79,'energy derivative',/,t3,'l',t8,
     +       'energy',t26,'value',t39,'derivative',t53,'nodes',t68,
     +       'value',t81,'derivative',t95,'nodes',t107,'norm',t119,
     +       'wronskian')
 8002 FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

      IF (film) THEN
         DO ivac = 1,2
            vz0(ivac) = vz(nmz,ivac)
         ENDDO
      ENDIF
c+q_sl
      IF ((ndir.EQ.-3).AND.dos) THEN
         IF (odi%d1) STOP 'layer-resolved feature does not work with 1D'
         CALL slab_dim(
     >                  natd,ntypd,
     >                  ntype,pos,neq,
     <                  nsld)
            ALLOCATE ( nmtsl(ntypd,nsld),nslat(natd,nsld) )
            ALLOCATE ( zsl(2,nsld),volsl(nsld) )
            ALLOCATE ( volintsl(nsld) )
            CALL slabgeom(
     >                    natd,ntypd,nsld,
     >                    ntype,pos,z1,neq,area,volmts,
     <                    nsl,zsl,nmtsl,nslat,volsl,volintsl)
c
        ALLOCATE ( qintsl(nsld,neigd,nkptd,jspd))
        ALLOCATE ( qmtsl(nsld,neigd,nkptd,jsp_start:jsp_end))
        ALLOCATE ( orbcomp(neigd,23,natd,nkptd,jsp_start:jsp_end) )
        ALLOCATE ( qmtp(neigd,natd,nkptd,jsp_start:jsp_end) )
        IF (.not.film) qvac(:,:,:,jspin) = 0.0
      ENDIF
c-q_sl
c
c   Open direct access file for cdninf-data if on cray T3E! This is now
c   always done, even if not on T3E no vacdos/dosinp file is written anymore
c     
      IF ((jspin.EQ.1).AND.(dos.or.vacdos.or.cdinf).
     +                 AND.(.not.semic)) THEN
c
c      calculate record length & open
c
         reclength= neigd*6 + 4*ntypd*neigd + neigd*layerd*2*jspd
         reclength= 8*(reclength + 4*nstars*neigd*layerd +5)
         IF (l_mcd) reclength= reclength + 8*(3*ncored*neigd*ntypd)
         IF (l_p0) THEN
           OPEN (84,file='tmp_dos',status='unknown',form='unformatted',
     +                                 access='direct',recl=reclength)
c+new_sl
           IF ((ndir.EQ.-3).AND.l_p0) THEN
             reclength2= 8*(5+neigd*(3+2*nsld+24*natd))
             OPEN (129,file='tmp_Ek',status='unknown',form='unformatted'
     +                                 ,access='direct',recl=reclength2)
           ENDIF
c-new_sl
         ENDIF
      ENDIF
c+dw  if on T3E assign to tmp_vacwave file
      IF (nstm.EQ.3) THEN
c
        reclength_vw= 8*( neigd + 4*n2max_nstm3*neigd + 4)
        OPEN (88,file='tmp_vacwave',status='unknown',form='unformatted',
     +     access='direct',recl=reclength_vw)
      ENDIF
c     open vacwave!
      if (nstm.EQ.3)  OPEN (87,file='vacwave',form='formatted',
     +     status='unknown')
c-dw
c+t3e
c
c-->   loop over k-points: each can be a separate task
c
      IF (nkpt < isize) THEN
        l_evp = .true.
        IF (l_mcd) THEN
           mcd(:,:,:) = 0.0
        ENDIF
        ener(:,:,:) = 0.0
        sqal(:,:,:) = 0.0
        qal(:,:,:,:) = 0.0
        enerlo(:,:,:) = 0.0
        sqlo(:,:,:) = 0.0
      ELSE
        l_evp = .false.
      ENDIF
      ALLOCATE ( we(neigd) )
      i_rec = 0 ; n_rank = 0
      DO 200 ikpt = 1,nkpt
        i_rec = i_rec + 1
         IF ((mod(i_rec-1,isize).EQ.irank).OR.l_evp) THEN
         skip_t = 0
         DO i = 1,skip_tt
           IF (mod(i-1,n_size).EQ.n_rank) skip_t = skip_t + 1
         ENDDO
         skip_t = skip_tt
c-t3e
         DO i = 1,neigd
            we(i) = 0.
         ENDDO
c--->    determine number of occupied bands and set weights (we)
         IF (semic) THEN
            noccbd = nint(zelec(nw)/2.0)
            IF (l_noco) noccbd = nint(zelec(nw))
            DO 70 i = 1,noccbd
               we(i) = wk
   70       CONTINUE
         ELSE
            noccbd = 0
            DO 80 i = 1,neigd ! nbands
               we(i) = w(n_bands(n_rank)+i,ikpt,jspin)
               IF (l_noco) we(i) = w(i,ikpt,1)
               IF ((we(i).GE.1.e-8).OR.pallst) THEN
                  noccbd = noccbd + 1
               ELSE
                  we(i)=0.0
               ENDIF 
   80       CONTINUE
         END IF
!
! -> Gu test: distribute ev's among the processors...
!
         IF (l_evp.AND.(isize.GT.1)) THEN
           IF (dos) THEN
             noccbd_l = CEILING( real(n_bands(1)) / isize )
             n_start = irank*noccbd_l + 1
             n_end   = min( (irank+1)*noccbd_l , n_bands(1) )
           ELSE
             noccbd_l = CEILING( real(noccbd) / isize )
             n_start = irank*noccbd_l + 1
             n_end   = min( (irank+1)*noccbd_l , noccbd )
           ENDIF
           noccbd = n_end - n_start + 1
           we(1:noccbd) = we(n_start:n_end)

           IF (n_start > skip_tt) THEN
              skip_t  = 0
           ENDIF
           IF (n_end <= skip_tt) THEN
              skip_t  = noccbd
           ENDIF
           IF ((n_start <= skip_tt).AND.(n_end > skip_tt)) THEN
              skip_t  = mod(skip_tt,noccbd)
           ENDIF
         ELSE
           n_start = 1
           IF (dos) THEN
             noccbd_l = n_bands(1)
             n_end    = n_bands(1)
             noccbd   = n_bands(1)
           ELSE
             noccbd_l = noccbd
             n_end    = noccbd
           ENDIF
         ENDIF
         ALLOCATE ( z(nbasfcn,neigd) )
         CALL cdn_read(
     >                 lmaxd,ntypd,nlod,neigd,nvd,jspd,irank,isize,
     >                 ikpt,jspin,nbasfcn,nlotot,l_ss,l_noco,
     >                 nrec,i_rec,ntapwf,noccbd,n_start,n_end,
     <                 nmat,nv,ello,evdu,epar,kveclo,
     <                 k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time)
         IF (l_evp.AND.(isize.GT.1)) THEN
           eig(1:noccbd) = eig(n_start:n_end)
         ENDIF
c
         IF (nstm.EQ.3.AND.film) THEN
           nat = 0
           DO n = 1,ntype
             nat = nat + neq(n)
           ENDDO
           CALL nstm3(
     >          nv2d,nop,natd,nvd,nmzd,k1d,k2d,k3d,n3d,ikpt,nv(jspin),
     >          jspins,nop2,mx1,mx2,n2d,jspin,nmz,nvac,nkpt,nw,tworkf,
     >          area,delz,z1,wk,k1(1,jspin),k2(1,jspin),jspd,nwdd,nat,
     >          ig2,ig,nstr2,mrot,rgphs,amat,bmat,evac(1,jspin),vz,vz0,
     >          taual,
     <          gvac1d,gvac2d)
         END IF

         IF (noccbd.EQ.0) GO TO 199
c
c--->    if slice, only a certain bands are taken into account
c--->    in order to do this the coresponding eigenvalues, eigenvectors
c--->    and weights have to be copied to the beginning of the arrays
c--->    eig, z and we and the number of occupied bands (noccbd) has to
c--->    changed
         IF (slice) THEN
            IF (l_p0) WRITE (16,FMT=*) 'NNNE',nnne
            IF (l_p0) WRITE (16,FMT=*) 'kk',kk
            nslibd = 0
            IF (pallst) THEN
               DO i = 1,nbands
                  we(i) = wk
               ENDDO
            END IF   
            IF (kk.EQ.0) THEN
               IF (l_p0) THEN
                 WRITE (16,FMT='(a)') 'ALL K-POINTS ARE TAKEN IN SLICE'
                 WRITE (16,FMT='(a,i2)') ' slice: k-point nr.',ikpt
               ENDIF
               DO i = 1,nbands
                  IF (eig(i).GE.e1s .AND. eig(i).LE.e2s) THEN
                     nslibd = nslibd + 1
                     eig(nslibd) = eig(i)
                     we(nslibd) = we(i)
                     DO  j = 1, nmat ! nv(jspin) + nlotot
                        z(j,nslibd) = z(j,i)
                     ENDDO
                  END IF
               ENDDO
               IF (l_p0)
     +         WRITE (16,'(a,i3)') ' eigenvalues in slice:',nslibd
            ELSEIF (kk.EQ.ikpt) THEN
               IF (l_p0) 
     +         WRITE (16,FMT='(a,i2)') ' slice: k-point nr.',ikpt
               IF ((e1s.EQ.0.0) .AND. (e2s.EQ.0.0)) THEN
                 IF (l_p0) 
     +           WRITE (16,FMT='(a,i5,f10.5)') 'slice: eigenvalue nr.',
     +              nnne,eig(nnne)
                  nslibd = nslibd + 1
                  eig(nslibd) = eig(nnne)
                  we(nslibd) = we(nnne)
                  DO j = 1, nmat ! nv(jspin) + nlotot
                     z(j,nslibd) = z(j,nnne)
                  ENDDO
               ELSE
                  DO i = 1,nbands
                     IF (eig(i).GE.e1s .AND. eig(i).LE.e2s) THEN
                        nslibd = nslibd + 1
                        eig(nslibd) = eig(i)
                        we(nslibd) = we(i)
                        DO j = 1, nmat ! nv(jspin) + nlotot
                          z(j,nslibd) = z(j,i)
                       ENDDO
                     END IF
                  ENDDO
                  IF (l_p0) 
     +            WRITE (16,FMT='(a,i3)')' eigenvalues in slice:',nslibd
               END IF
            END IF
            noccbd = nslibd
            IF (nslibd.EQ.0) GO TO 199 ! 200
         END IF ! slice

c--->    in normal iterations the charge density of the unoccupied
c--->    does not need to be calculated (in pwden, vacden and abcof)
         IF (dos.AND. .NOT.(l_evp.AND.(isize.GT.1)) ) THEN
             noccbd=nbands
         ENDIF
        

c     ----> add in spin-doubling factor
         DO 90 i = 1,noccbd
            we(i) = 2.*we(i)/jspins
   90    CONTINUE

c---> pk non-collinear
c--->    valence density in the interstitial and vacuum region
c--->    has to be called only once (if jspin=1) in the non-collinear 
c--->    case
c     ----> valence density in the interstitial region
         IF (.NOT.((jspin.EQ.2) .AND. l_noco)) THEN
         CALL cpu_time(time1)
         CALL pwden(
     >        nvd,nbasfcn,neigd,n3d,nkptd,jspd,noccbd,k1d,k2d,k3d,
     >        kq1d,kq2d,kq3d,mx1,mx2,mx3,ig,rgphs,jspins,isize,ikpt,
     >        jspin,l_noco,l_ss,k1,k2,k3,noccbd,nv,nstr,kq1_fft,kq2_fft,
     >        kq3_fft,nq3_fft,kmxq_fft,igq_fft,igfft,pgfft,omtil,we,z,
     >        ntypd,natd,nop,ntype,z1,vol,volint,l_f,nlotot,invtab,
     >        symor,cdinf,dos,vacdos,kv3,mrot,neq,ig,ig2,bbmat,
     >        taual,volmts,tau,bmat,rmt,sk3,eig,bkpt,lmax,pos,lmaxd,odi,
     +        qpw,cdom,qis,force,f_b8)
         CALL cpu_time(time2)
         cp_time(2) = cp_time(2) + time2 - time1
         ENDIF
c+new
c--->    charge of each valence state in this k-point of the SBZ
c--->    in the layer interstitial region of the film
c
         IF (dos.AND.(ndir.EQ.-3))  THEN
            IF (.NOT.((jspin.EQ.2) .AND. l_noco)) THEN
              CALL q_int_sl(
     >                    nvd,nbasfcn,neigd,n3d,natd,nop,invtab,
     >                    k1d,k2d,k3d,ntypd,ntype,neq,nq3,
     >                    volmts,volsl,volintsl,taual,symor,tau,mrot,
     >                    omtil,mx1,mx2,mx3,nstr,rgphs,ig,ig2,
     >                    z,noccbd,nv(jspin),k1(1,jspin),k2(1,jspin),
     >                    k3(1,jspin),rmt,sk3,bmat,kv3,
     >                    nsl,zsl,nmtsl,odi,
     <                    qintsl(1,1,ikpt,jspin))
c
            ENDIF
           ENDIF
c-new c
c--->    valence density in the vacuum region
         IF (film) THEN
            IF (.NOT.((jspin.EQ.2) .AND. l_noco)) THEN
            CALL cpu_time(time1)
            CALL vacden(
     >           nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,nvd,jspd,
     >           nkptd,neigd,layerd,nvac,jspins,nmz,nmzxy,
     >           mx1,mx2,z1,area,delz,omtil,bmat,ikpt,wk,
     >           nstr2,ig2,rgphs,ig,gvac1d,gvac2d,nlotot,
     >           we,ikpt,jspin,l_noco,l_ss,npotmatfile,vz,vz0,
     >           noccbd,z,bbmat,bkpt,k1,k2,k3,nv,nbasfcn,nmat,
     >           layers,izlay,vacdos,integ,evac,eig,
     >           tworkf,nstars,nstm,starcoeff,locx,locy,qss,
     >           odi,ods,sk2,phi2,
     X           rhtxy,rht,qvac,qvlay,emax_nstm3,n2max_nstm3,
     <           qstars,cdomvz,cdomvxy)
            CALL cpu_time(time2)
            cp_time(3) = cp_time(3) + time2 - time1
            ENDIF
c--->       perform Brillouin zone integration and summation over the
c--->       bands in order to determine the vacuum energy parameters.
            DO ispin = jsp_start,jsp_end
               DO ivac = 1,nvac
                  DO ie = 1,noccbd
                     pvac(ivac,ispin) = pvac(ivac,ispin) + eig(ie) *
     +                                qvac(ie,ivac,ikpt,ispin)*we(ie)
                     svac(ivac,ispin) = svac(ivac,ispin)
     +                              + qvac(ie,ivac,ikpt,ispin)*we(ie)
                  ENDDO
               ENDDO
            ENDDO
         END IF

c--->    valence density in the atomic spheres
c--->    construct a(tilta) and b(tilta)
         IF (l_mperp) THEN
            ALLOCATE ( acof(noccbd,0:lmd,natd,jspd),
! Deallocated before call to sympsi
     +                 bcof(noccbd,0:lmd,natd,jspd),                
     +                 ccof(-llod:llod,noccbd,nlod,natd,jspd) )
         ELSE
            ALLOCATE ( acof(noccbd,0:lmd,natd,jspin:jspin),
     +                 bcof(noccbd,0:lmd,natd,jspin:jspin),
     +                 ccof(-llod:llod,noccbd,nlod,natd,jspin:jspin) )
         ENDIF

         DO ispin = jsp_start,jsp_end
            CALL cpu_time(time1)
            IF (l_f) THEN
c+for
               ALLOCATE (e1cof(noccbd,0:lmaxd*(lmaxd+2),natd),          ! Deallocated after call to force_a21
     +              e2cof(noccbd,0:lmaxd*(lmaxd+2),natd),
     +              acoflo(-llod:llod,noccbd,nlod,natd),
     +              bcoflo(-llod:llod,noccbd,nlod,natd),
     +              aveccof(3,noccbd,0:lmaxd*(lmaxd+2),natd),
     +              bveccof(3,noccbd,0:lmaxd*(lmaxd+2),natd),
     +              cveccof(3,-llod:llod,noccbd,nlod,natd) )

               CALL to_pulay(
     >           lmaxd,ntypd,neigd,noccbd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nmat,nlotot,invtab,
     >           l_noco,l_ss,lmax,ntype,nv,omtil,bkpt,
     >           k1,k2,k3,neq,rmt,mrot,ngopr,invsat,invsatnr,
     >           taual,bmat,bbmat,z,noccbd,eig,us(0,1,ispin),
     >           dus(0,1,ispin),uds(0,1,ispin),duds(0,1,ispin),
     >           ddn(0,1,ispin),ulos(1,1,ispin),uulon(1,1,ispin),
     >           dulon(1,1,ispin),dulos(1,1,ispin),llo,nlo,
     >           l_dulo,lapw_l,kveclo,l_geo,ispin,alph,beta,qss,odi,ods,
     <           acof(1,0,1,ispin),bcof(1,0,1,ispin),
     <           e1cof,e2cof,aveccof,bveccof,
     <           ccof(-llod,1,1,1,ispin),acoflo,bcoflo,cveccof)
c-for
            ELSE
               CALL abcof(
     >           lmaxd,ntypd,neigd,noccbd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngopr,taual,neq,lmax,rmt,omtil,
     >           bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,noccbd,z,
     >           us(0,1,ispin),dus(0,1,ispin),uds(0,1,ispin),
     >           duds(0,1,ispin),ddn(0,1,ispin),invsat,invsatnr,
     >           ulos(1,1,ispin),uulon(1,1,ispin),dulon(1,1,ispin),
     >           dulos(1,1,ispin),llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,ispin,alph,beta,qss,kveclo,odi,ods,
     <           acof(1,0,1,ispin),bcof(1,0,1,ispin),
     <           ccof(-llod,1,1,1,ispin))


            ENDIF
            CALL cpu_time(time2)
            cp_time(4) = cp_time(4) + time2 - time1

            IF ( n_u.GT.0 ) THEN
            CALL n_mat(
     >                 lmaxd,ntypd,neigd,noccbd,natd,nop,llod,nlod,
     >                 ngopr,neq,lmax,noccbd,ddn(0,1,ispin),we,ntype,
     >                 lda_u,n_u,uulon(1,1,ispin),dulon(1,1,ispin),
     >                 llo,nlo,uloulopn(1,1,1,ispin),invsat,invtab,lmd,
     >                 d_wgn,acof(1,0,1,ispin),bcof(1,0,1,ispin),
     >                 ccof(-llod,1,1,1,ispin),invarind,invarop,
     <                 n_mmp)
            ENDIF
c
c--->       perform Brillouin zone integration and summation over the
c--->       bands in order to determine the energy parameters for each
c--->       atom and angular momentum
c
            IF (.not.slice) THEN
              CALL eparas(
     >        llod,noccbd,nlod,natd,neigd,ntypd,lmaxd,lmd,isize,
     >        ikpt,ntype,neq,noccbd,we,eig,ccof(-llod,1,1,1,ispin),
     >        nlo,llo,invsat,skip_t,l_evp,
     >        acof(1,0,1,ispin),bcof(1,0,1,ispin),ddn(0,1,ispin),
     >        ncore,nstd,l_mcd,m_mcd,uulon(1,1,ispin),
     >        dulon(1,1,ispin),uloulopn(1,1,1,ispin),
     <        enerlo(1,1,ispin),sqlo(1,1,ispin),
     <        ener(0,1,ispin),sqal(0,1,ispin),
     <        qal(0,1,1,ispin),mcd)

              IF ( l_mperp .AND. (ispin == jsp_end) ) THEN
                CALL qal_21(
     >                      llod,nlod,natd,neigd,ntypd,lmaxd,lmd,
     >                      jspins,ntype,neq,noccbd,we,ccof,nlo,llo,
     >                     alph,beta,acof,bcof,mt21,lo21,uloulopn21,
     <                     qal,qmat)
              ENDIF
            ENDIF
c
c+new
c--->    layer charge of each valence state in this k-point of the SBZ
c--->    from the mt-sphere region of the film
c
             IF (dos.AND.(ndir.EQ.-3))  THEN
               write(*,*) 'q_mt_sl'
               CALL q_mt_sl(
     >             llod,noccbd,nlod,natd,neigd,ntypd,lmaxd,lmd,nsld,
     >             ikpt,ntype,neq,noccbd,ccof(-llod,1,1,1,ispin),
     >             nlo,llo,invsat,skip_t,noccbd,
     >             acof(1,0,1,ispin),bcof(1,0,1,ispin),ddn(0,1,ispin),
     >             uulon(1,1,ispin),
     >             dulon(1,1,ispin),uloulopn(1,1,1,ispin),
     >             nmtsl,nsl,lmax,
     <             qmtsl(1,1,ikpt,ispin))
c

               INQUIRE (file='orbcomprot',exist=l_orbcomprot)
               IF (l_orbcomprot) THEN                           ! rotate ab-coeffs
                  CALL abcrot2(
     >                         ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,
     >                         ntype,neq,noccbd,lmax,nlo,llo,
     X                         acof(1,0,1,ispin),bcof(1,0,1,ispin),
     X                         ccof(-llod,1,1,1,ispin))
               ENDIF 

               CALL orb_comp(
     >                  noccbd,lmd,lmaxd,natd,ntypd,nlod,llod,
     >                  ntype,noccbd,neq,ddn(0,1,ispin),
     >                  uulon(1,1,ispin),dulon(1,1,ispin),
     >                  uloulopn(1,1,1,ispin),
     >                  acof(1,0,1,ispin),bcof(1,0,1,ispin),
     >                  ccof(-llod,1,1,1,ispin),nlo,llo,neigd,
     <                  orbcomp(1,1,1,ikpt,ispin),
     <                  qmtp(1,1,ikpt,ispin))
c
             ENDIF
c-new
c--->          set up coefficients for the spherical and
               CALL cpu_time(time1)
               CALL rhomt(
     >                 lmaxd,noccbd,natd,ntypd,we,noccbd,ntype,neq,lmax,
     >                 acof(1,0,1,ispin),bcof(1,0,1,ispin),
     X                 uu(0,1,ispin),dd(0,1,ispin),du(0,1,ispin))
               CALL cpu_time(time2)
               cp_time(6) = cp_time(6) + time2 - time1
c+soc            
               IF (l_soc) THEN
                  CALL orbmom(
     >                ntypd,lmaxd,natd,noccbd,lmd,ntype,neq,noccbd,lmax,
     >                we,acof(1,0,1,ispin),bcof(1,0,1,ispin),
     >                nlod,llod,nlo,llo,ccof(-llod,1,1,1,ispin),
     X                orb(0,-lmaxd,1,ispin),orbl(1,-llod,1,ispin),
     X                orblo(1,1,-llod,1,ispin) )
               ENDIF
c     -soc
c--->          non-spherical m.t. density
               CALL cpu_time(time1)
               CALL rhonmt(
     >           lmaxd,ntypd,memd,nlhd,ntypsd,natd,noccbd,
     >           we,noccbd,ntype,nsymt,neq,lmax,
     >           clnu,nmem,mlh,nlh,llh,ntypsy,
     >           acof(1,0,1,ispin),bcof(1,0,1,ispin),
     X           uunmt(0,1,1,ispin),ddnmt(0,1,1,ispin),
     X           udnmt(0,1,1,ispin),dunmt(0,1,1,ispin))
               CALL cpu_time(time2)
               cp_time(7) = cp_time(7) + time2 - time1

c--->          set up coefficients of the local orbitals and the 
c--->          flapw - lo cross terms for the spherical and 
c--->          non-spherical mt density
               CALL rhomtlo(
     >           noccbd,natd,lmd,nlod,llod,ntypd,ntype,
     >           neq,noccbd,we,acof(1,0,1,ispin),bcof(1,0,1,ispin),
     >           ccof(-llod,1,1,1,ispin),nlo,llo,invsat,
     X           aclo(1,1,ispin),bclo(1,1,ispin),cclo(1,1,1,ispin))
c
               CALL rhonmtlo(
     >           natd,lmd,nlod,memd,noccbd,ntypd,nlhd,lmaxd,
     >           ntype,neq,lmax,noccbd,we,acof(1,0,1,ispin),
     >           bcof(1,0,1,ispin),ccof(-llod,1,1,1,ispin),llod,
     >           clnu,nmem,mlh,nlh,llh,ntypsy,ntypsd,nlo,llo,
     X           acnmt(0,1,1,1,ispin),bcnmt(0,1,1,1,ispin),
     X           ccnmt(1,1,1,1,ispin))
c
               IF (l_f) THEN
c+for
#ifndef CPP_APW
                  CALL force_a12(
     >              lmaxd,ntypd,noccbd,natd,nop,
     >              jspd,mrot,ngopr,lmax,invarop,invarind,
     >              ntype,neq,rmt,invtab,multab,amat,bmat,
     >              we,ispin,noccbd,us,uds,acof(1,0,1,ispin),
     >              bcof(1,0,1,ispin),e1cof,e2cof,
     >              lmd,nlod,llod,nlo,llo,acoflo,bcoflo,l_geo,
     >              odi,ods,
     X              force,f_a12)
#endif
                  CALL force_a21(
     >              lmaxd,ntypd,neigd,noccbd,natd,nop,nlod,llod,loplod,
     >              lmd,jspd,mrot,ngopr,lmax,invarop,invarind,ntype,neq,
     >              invtab,multab,amat,bmat,we,ispin,nw,epar(0,1,ispin),
     >              noccbd,eig,ddn(0,1,ispin),l_geo,lda_u,
     >              nlo,llo,nlol,lo1l,dulon(1,1,ispin),uulon(1,1,ispin),
     >              uloulopn(1,1,1,ispin),acof(1,0,1,ispin),
     >              bcof(1,0,1,ispin),ccof(-llod,1,1,1,ispin),
     >              aveccof,bveccof,cveccof,
     >              us(0,1,ispin),uds(0,1,ispin),dus(0,1,ispin),
     >              duds(0,1,ispin),ulos(1,1,ispin),dulos(1,1,ispin),
     >              rmt,odi,ods,
     X              force,f_a21,f_b4)
c     -for
              
            ENDIF
            IF (l_f) THEN
               DEALLOCATE (e1cof,e2cof,aveccof,bveccof)
               DEALLOCATE (acoflo,bcoflo,cveccof)
            ENDIF
c--->    end loop over ispin
         ENDDO
         
         IF (l_mperp) THEN
            CALL rhomt21(
     >                   lmaxd,noccbd,natd,ntypd,jspd,
     >                   we,noccbd,ntype,neq,lmax,acof,bcof,
     >                   nlod,llod,nlo,llo,ccof,
     X                   mt21,lo21,uloulop21)
            IF (l_fmpl) THEN
               CALL rhonmt21(
     >                       lmaxd,llpd,ntypd,memd,nlhd,ntypsd,natd,
     >                       noccbd,jspd,we,noccbd,ntype,nsymt,neq,lmax,
     >                       clnu,nmem,mlh,nlh,llh,ntypsy,
     >                       acof,bcof,
     X                       uunmt21,ddnmt21,udnmt21,dunmt21)
            ENDIF
         ENDIF

         DEALLOCATE (acof,bcof,ccof)
c
  199    CONTINUE
         IF ((dos .OR. vacdos .OR. cdinf) .AND.(.not.semic) ) THEN
            CALL cpu_time(time1)
c
c--->    calculate charge distribution of each state (l-character ...)
c--->    and write the information to the files dosinp and vacdos
c--->    for dos and bandstructure plots
c

c--dw    parallel writing of vacdos,dosinp....
c        write data to direct access file first, write to formated file later by PE 0 only!
c--dw    since z is no longer an argument of cdninf sympsi has to be called here!
c
            CALL cotra3(bkpt,cartk,bmat)
            IF (ndir.GT.0) THEN
               CALL sympsi(
     >                     bkpt,nv(jspin),k1(1,jspin),k2(1,jspin),
     >                     k3(1,jspin),nop,neigd,nbands,nvd,mrot,amat,
     >                     nbasfcn,z,eig,l_noco,l_soc,
     <                     ksym,jsym)
            END IF
c
c--dw   now write k-point data to tmp_dos
c
            IF (l_p0) THEN
              IF ( .not.l_mcd ) THEN
                DO ispin = jsp_start,jsp_end
                  WRITE (84,rec=nkpt*(ispin-1)+ikpt) bkpt,wk,nbands,eig,
     +                            qal(:,:,:,ispin),qvac(:,:,ikpt,ispin),
     +                        qis(:,ikpt,ispin),qvlay(:,:,:,ikpt,ispin),
     +                                                  qstars,ksym,jsym
                ENDDO
              ELSE
                WRITE (84,rec=nkpt*(jspin-1)+ikpt) bkpt,wk,nbands,eig,
     +                          qal(:,:,:,jspin),qvac(:,:,ikpt,jspin),
     +                      qis(:,ikpt,jspin),qvlay(:,:,:,ikpt,jspin),
     +                              qstars,ksym,jsym,mcd(:,1:ncored,:)
              ENDIF ! not.l_mcd
              CALL cpu_time(time2)
              cp_time(5) = cp_time(5) + time2 - time1
c+new_sl
              IF (ndir.EQ.-3) THEN
                IF ( .not.l_mcd ) THEN
                   DO ispin = jsp_start,jsp_end
                      WRITE (129,rec=nkpt*(ispin-1)+ikpt) bkpt,wk,
     +                            nbands,eig,qvac(:,:,ikpt,ispin),
     +                            qintsl(:,:,ikpt,ispin),
     +                            qmtsl(:,:,ikpt,ispin),
     +                            orbcomp(:,:,:,ikpt,ispin),
     +                            qmtp(:,:,ikpt,ispin)
                   ENDDO
                ELSE
                   WRITE (6,FMT=500)
                   WRITE (16,FMT=500)
 500               FORMAT('A calculation of an orbital composition
     +                   of film states  (s,px,py,pz,....)
     +                   for the l_mcd  case is not implemented.')
                ENDIF ! not.l_mcd
              ENDIF ! ndir == -3
            ENDIF ! l_p0
#ifdef CPP_MPI
            DO ispin = jsp_start,jsp_end
              CALL mpi_col_dos(
     >                       l_evp,noccbd,noccbd_l,
     >                       irank,isize,jspd,neigd,lmaxd,ntypd,
     >                       nsld,natd,ncored,nstd,
     >                       nstars,layerd,l_mcd,ndir,nkpt,ispin,ikpt,
     >                       bkpt,wk,nbands,eig,qal(0,1,1,ispin),
     >                       qvac(1,1,ikpt,ispin),qis(1,nkpt,ispin),
     >                       qvlay(:,:,:,ikpt,ispin),qstars,ksym,
     >                       jsym,m_mcd,qintsl(:,:,ikpt,ispin),
     >                       qmtsl(:,:,ikpt,ispin),qmtp(:,:,ikpt,ispin),
     >                       orbcomp(:,:,:,ikpt,ispin))
            ENDDO
#endif
c-new_sl
         ENDIF

c--->  end of loop over PE's
         DEALLOCATE (z)
       ENDIF 
c---> end of k-point loop
  200 CONTINUE
      nrec = nrec + nkpt*n_size
      DEALLOCATE (we,f,g,us,dus,duds,uds,ddn)
c+t3e
#ifdef CPP_MPI
      CALL cpu_time(time1)
      DO ispin = jsp_start,jsp_end
          CALL mpi_col_den(
     >          irank,isize,ntypsd,ntype,nlod,natd,llod,odi,
     >          lmaxd,ntypd,nlhd,n3d,n2d,nmzd,nmzxyd,
     >          jspins,nq2,nq3,nmz,nmzxy,lmax,nlh,ntypsy,
     >          film,l_f,l_soc,l_mperp,l_fmpl,ispin,nw,nlo,neq,n_u,llpd,
     X          rhtxy(1,1,1,ispin),rht(1,1,ispin),qpw(1,ispin),
     X          ener(0,1,ispin),sqal(0,1,ispin),force(1,1,ispin),
     X          svac(1,ispin),pvac(1,ispin),uu(0,1,ispin),dd(0,1,ispin),
     X          du(0,1,ispin),uunmt(0,1,1,ispin),ddnmt(0,1,1,ispin),
     X          udnmt(0,1,1,ispin),dunmt(0,1,1,ispin),sqlo(1,1,ispin),
     X          aclo(1,1,ispin),bclo(1,1,ispin),cclo(1,1,1,ispin),
     X          acnmt(0,1,1,1,ispin),bcnmt(0,1,1,1,ispin),
     X          ccnmt(1,1,1,1,ispin),enerlo(1,1,ispin),
     X          orb(0,-lmaxd,1,ispin),orbl(1,-llod,1,ispin),
     X          orblo(1,1,-llod,1,ispin),mt21,lo21,uloulop21,
     X          uunmt21,ddnmt21,udnmt21,dunmt21,
     X          cdom,cdomvz,cdomvxy,l_noco,n_mmp)
      ENDDO
      CALL cpu_time(time2)
      cp_time(9) = cp_time(9) + time2 - time1

c     all PE close 88  -> cache is written to disk!
      IF (nstm.EQ.3) CLOSE(88)
#endif 
c--dw now close tmp_dos and open it again to ensure that pe0 gets data from all pe's
c
       IF (((jspin.eq.jspins).OR.l_mperp) .AND. (dos.or.vacdos.or.cdinf) 
     +                                         .AND. (.NOT.semic)) THEN
        CLOSE (84)

         IF ((ndir.EQ.-3).AND.l_p0) CLOSE (129)
         reclength= neigd*6 + 4*ntypd*neigd + neigd*layerd*2*jspd
         reclength= 8*(reclength + 4*nstars*neigd*layerd +5)
         IF (l_mcd) reclength= reclength + 8*(3*ncored*neigd*ntypd)
         OPEN (84,file='tmp_dos',status='old',form='unformatted',
     +         access='direct',recl=reclength)
c+new_sl
         IF ((ndir.EQ.-3).AND.l_p0) THEN
           reclength2= 8*(5+neigd*(3+2*nsld+24*natd))
           OPEN (129,file='tmp_Ek',status='old',form='unformatted',
     +                             access='direct',recl=reclength2)
         ENDIF
c
         IF (l_p0) THEN
            CALL doswrite(
     >                    neigd,nkptd,nvd,ntypd,jspd,layerd,
     >                    jspins,semic,dos,vacdos,cdinf,film,
     >                    nkpt,ntype,layers,nvac,slice,l_noco,nop2,
     >                    ndir,invs,zrfs,bmat,neq,izlay,
     >                    n2max_nstm3,reclength_vw,nstm,tworkf,
     >                    nstars,starcoeff,l_mcd,ncored,ncore,e_mcd,
     >                    e1_dos,e2_dos,efermi,sig_dos,nsld,natd,odi)
c+new_sl
            IF (dos.AND.(ndir.EQ.-3))  THEN
                CALL Ek_write_sl(
     >                           neigd,nkptd,ntypd,jspd,layerd,natd,
     >                           nsld,nwdd,jspins,nkpt,ntype,nvac,jspin,
     >                           nw,invs,zrfs,bmat,neq,
     >                           nsl,nslat)
            ENDIF
c-new_sl
         ENDIF
      ENDIF

      IF (l_p0) THEN
         CALL cdnmt(
     >        jspd,ntypd,natd,lmaxd,ntypsd,llpd,nlhd,jmtd,nlod,llod,
     >        l_soc,l_mperp,l_fmpl,ntype,jsp_start,jsp_end,sfp,
     >        lmax,nlh,llh,neq,ntypsy,jri,nlo,llo,l_dulo,ulo_der,
     >        dx,rmsh,epar,ello,vr,uu,du,dd,uunmt,udnmt,dunmt,ddnmt,
     >        uulon,dulon,uloulopn,aclo,bclo,cclo,acnmt,bcnmt,ccnmt,
     >        orb,orbl,orblo,mt21,lo21,uloulopn21,uloulop21,
     >        uunmt21,ddnmt21,udnmt21,dunmt21,
     <        chmom,clmom,
     X        qa21,cp_time(8),rho)

         DO ispin = jsp_start,jsp_end
           WRITE (6,*) 'Energy Parameters for spin:',ispin
           IF (.not.slice) THEN
           CALL mix_enpara(
     >                     ntype,nlod,lmaxd,nmzd,jmtd,lepr,film,
     >                     nvac,enmix(ispin,nw),nlo,llo,jri,dx,rmt,
     >                     vr(1,1,ispin),vz,pvac(1,ispin),
     >                     svac(1,ispin),lchg_v(1,ispin,nw),
     >                     ener(0,1,ispin),sqal(0,1,ispin),
     >                     lchange(0,1,ispin,nw),enerlo(1,1,ispin),
     >                     sqlo(1,1,ispin),llochg,l_dulo,
     X                     ello0(1,1,ispin),el0(0,1,ispin,nw),
     X                     evac0(1,ispin,nw))
           CALL w_enpara(
     >                   lmaxd,nlod,ntype,nw,jspin,film,nlo,
     >                   skiplo,ello0(1,1,ispin),el0(0,1,ispin,nw),
     >                   evac0(1,ispin,nw),lchange(0,1,ispin,nw),
     >                   llochg,lchg_v(1,ispin,nw),enmix(ispin,nw),16)
           ENDIF

c--->      check continuity of charge density
           IF (cdinf .AND. .NOT.semic) THEN
               WRITE (6,FMT=8210) ispin
 8210          FORMAT (/,5x,'check continuity of cdn for spin=',i2)
               IF (film .AND. .NOT.odi%d1) THEN
c--->             vacuum boundaries
                  npd = min(nspd,25)
                  CALL points(xp,npd)
                  DO ivac = 1,nvac
                     sign = 3. - 2.*ivac
                     DO j = 1,npd
                        xp(3,j) = sign*z1/amat(3,3)
                     ENDDO
                     CALL checkdop(
     >                    xp,npd,0,0,ivac,1,ispin,.true.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    qpw,rho,rhtxy,rht,odi,ods)
                  ENDDO
               ELSE IF (odi%d1) THEN
c-odim
                  npd = min(nspd,25)
                  CALL cylpts(xp,npd,z1)
                     CALL checkdop(
     >                    xp,npd,0,0,ivac,1,ispin,.true.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    qpw,rho,rhtxy,rht,odi,ods)
c+odim
               END IF
c--->          m.t. boundaries
               nat = 1
               DO n = 1,ntype
                  CALL sphpts(xp,nspd,rmt(n),pos(1,nat))
                  CALL checkdop(
     >                    xp,nspd,n,nat,0,-1,ispin,.true.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    qpw,rho,rhtxy,rht,odi,ods)
                  nat = nat + neq(n)
               ENDDO
           ENDIF
c+for
c--->      forces of equ. A8 of Yu et al.
           IF ((l_f).AND.(.NOT.semic)) THEN
              CALL force_a8(
     >                 jmtd,memd,nlhd,ntypd,jspd,ntypsd,lmaxd,l_geo,
     >                 jri,mlh,nlh,nmem,ntype,ispin,ntypsy,natd,
     >                 dx,rmsh,clnu,llh,vrfull,rho,neq,
     >                 f_a12,f_a21,f_b4,f_b8,
     X                 force)
           ENDIF
c-for

         ENDDO ! end of loop ispin = jsp_start,jsp_end
      ENDIF ! end of (l_p0)
c+t3e
      IF (ndir.EQ.-3)  DEALLOCATE  (qintsl,qmtsl,orbcomp,qmtp)
      IF (ndir.EQ.-3)  DEALLOCATE  (nmtsl,nslat,zsl,volsl,volintsl)
      DEALLOCATE (orb,orbl,orblo)
      IF (l_f)   DEALLOCATE (f_a12,f_a21,f_b4)
      DEALLOCATE (f_b8)
      DEALLOCATE ( uunmt21,udnmt21,dunmt21,ddnmt21 )
      DEALLOCATE ( mt21,lo21,uloulop21,uloulopn21,qmat )
      DEALLOCATE ( m_mcd,mcd )
      DEALLOCATE ( ulos,dulos,uulon,dulon,uloulopn,kveclo )
      DEALLOCATE ( uu,dd,du,uunmt,ddnmt,dunmt,udnmt,nv,enerlo,sqlo )
      DEALLOCATE ( aclo,bclo,cclo,acnmt,bcnmt,ccnmt,gvac1d,gvac2d )
      DEALLOCATE ( k1,k2,k3,jsym,ksym,qal,sqal,ener,svac,pvac,qstars )

      IF ((jsp_end.EQ.jspins).AND.(.NOT.semic)) THEN
        IF ((ndir.EQ.-3).AND.l_p0) CLOSE (129)
#ifdef CPP_MPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
        IF ((dos.or.vacdos).and.(ndir.ne.-2)) STOP 'DOS OK'
        IF (nstm.EQ.3) STOP 'VACWAVE OK'
      ENDIF
      END SUBROUTINE cdnval
      END MODULE m_cdnval
