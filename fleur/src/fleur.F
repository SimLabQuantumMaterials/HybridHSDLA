      PROGRAM fleur
c     ***************************************************************
c
c     based on flapw7 (c.l.fu, m.weinert, e.wimmer):
c     full potential linearized augmented plane wave method for thin
c     films and superlattices (version 7 ---- general symmetry)
c     symmetry part       ---  e.wimmer
c     potential generator ---  c.l.fu,r.podloucky
c     matrix elements     ---  m.weinert
c     charge density      ---  c.l.fu
c                                c.l.fu        1987
c     2nd variation diagon.  --- r.-q. wu      1992
c     forces a la Yu et al   --- r.podloucky   1995
c     full relativistic core --- a.shick       1996
c     broyden mixing         --- r.pentcheva   1996
c     gga (pw91, pbe)        --- t.asada       1997
c     local orbitals         --- p.kurz        1997
c     automatic symmetry     --- w.hofer       1997
c     core tails & start     --- r.abt         1998
c     spin orbit coupling    --- a.shick,x.nie 1998
c     non-colinear magnet.   --- p.kurz        1999
c     one-dimensional        --- y.mokrousov   2002 
c     exchange parameters    --- m.lezaic      2004 
c     
c                       g.bihlmayer, s.bluegel 1999
c     ***************************************************************
c----------------------------------------
c this routine is the main PROGRAM
c
c flapw --+-- dimens      ... dimensioning routine(s)
c         +-- inped - rw_inp
c         +-- setup       ... readin & setup
c         +-- optional    ... plots, startdensity etc.
c        +-+
c        | +-- vgen       ... construct potential
c        | +-- eigen      ... calculate eigenvalues & vectors
c        |(+-- eigenso)       ... spin orbit coupling
c        | +-- fermie     ... calculate fermi energy
c        |(+-- force0)        ... initialize force
c        | +-- cdngen     ... construct charge density
c        |(+-- vgen)          ... output potential
c        |(+-- potdis)        ... potential distance
c        | +-- totale     ... calculate total energy
c        | +-- mix        ... mix densities
c        +-+
c----------------------------------------
c
c I/O file list 
c        2... time.info ...
c        5... inp       ...
c        6... out       ...
c        8... pottot    ...
c        9... nrp       ...
c       11... potcoul   ...
c       14... wkf2      ...
c       15... QGpsi,CLASS,NLAindx,SYMOPS,LATTC,gwa,abcoeff,ecore
c       16... inf       ...
c       17... cdnc      ...
c       18... plotin    ...
c       19... cdnout    ...
c       20... cdn_slice ...
c       22... dirofmag  ...
c       24... nocoinf   ...
c       25... potmat    ...
c       26... rhomat_inp & rhomat_out 
c       28... tmat      ...
c       29... tlst      ...
c       38... tmas      ...
c       40... enpara    ...
c       41... kpts      ...
c       43... force.dat ...
c       45... inp_new   ...
c       51... stars     ...
c       57... broyd     ...
c       59... broyd.x   ...
c       66... eig       ...
c       67... eig.soc   ...
c       71... cdn1      ...
c       79... cdn9      ...
c       85... dosinp    ...
c       86... vacdos    ...
c       87... tmp_vacdos ...
c       88... vacwave   ...
c       92... CBS_eigenvalues.1/2
c       93... CBS_states
c       99... qfix      ...
c      101... Smat...
c      102... Tmat..
c      110... scratch unit for GW
c      113... qpts      ...
c      114... jenerg    ...
c      115... jconst    ...
c      116... qptsinfo  ...
c      117... shells    ...
c      118... MCinp     ...
c      300... matelsso  (for ss+soc) 
c      301+XXX matelXXX (for ss+soc)
c
c----------------------------------------
#ifdef CPP_HPM
#include "/usr/local/hpm/LOCALinclude/f_hpm.h"
#endif
      USE m_types,     ONLY  : t_utype
      USE m_od_types 
      USE m_dimens
      USE m_inped
      USE m_setup
      USE m_optional
      USE m_vgen 
      USE m_rhodirgen
      USE m_vmatgen
      USE m_pldngen
      USE m_eigen
      USE m_eigenso
      USE m_fermie
      USE m_force0
      USE m_cdngen
      USE m_totale
      USE m_potdis
      USE m_mix
      USE m_cputime
      USE m_outtime
      USE m_jcoff
      USE m_jcoff2
#ifdef CPP_WANN
      USE m_wannier
#endif
#ifdef CPP_MPI
      USE m_mpi_bc_all,  ONLY : mpi_bc_all
#endif
#ifdef CPP_HDF
      USE m_hdf_tools,   ONLY : hdf_init, hdf_close
      USE m_eig66_hdf,   ONLY : openeig, closeeig
#endif
#ifdef DUMP_DATA
! need binmat module
      USE binmat
#endif
      IMPLICIT NONE
#ifdef DUMP_DATA_OLD
! variables for dumping stuff
      character(len=32) :: misc_name
#endif
c
c
c     scalar bmu changed to array bmu(ntypd) to enable calculation of
c     both ferromagnetic and antiferromagnetic CASE - r.pentcheva Jan'96
c     another possibility for the antiferromagnetic case ist to start
c     from the ferromagnetic charge density and flip the magnetic moments
c     of certain atoms: logical variable lflip and integer nflip(ntype)
c                                 r.pentcheva,kfa,Feb'96
C     .. Parameters ..
      INTEGER :: nop,k1d,k2d,k3d,n3d,n2d
      INTEGER :: kq1d,kq2d,kq3d
      INTEGER :: kxc1d,kxc2d,kxc3d
      INTEGER :: ntypd,natd,jmtd
      INTEGER :: ntypsd,nlhd,memd
      INTEGER :: lmaxd
      INTEGER :: jspd,nvacd
      INTEGER :: nmzd,nmzxyd
      INTEGER :: nspd
      INTEGER :: nvd,nv2d,nwdd,nkptd
      INTEGER :: neigd
      INTEGER :: msh,nstd
      INTEGER :: ncvd
      INTEGER :: layerd
      INTEGER :: nn2d,nn3d
      INTEGER :: nlod,llod
      INTEGER :: nqptd
C     ..
C     .. Local Scalars ..
      INTEGER icorr,izmax,izmin,jspins,krla,
     +        mx1,mx2,mx3,ndir,ndir_tmp,ng2,d1,d10,d100,aoff,
     +        ng3,ngz,nk1,nk2,nk3,nmz,nmzxy,nop2,nq2,nq3,nsymt,
     +        ntype,nvac,nbasfcn,TYP_REAL,TYP_COMPLEX
      INTEGER i,isec1,it,itmax,jspin,kcrel,lpr,n,n_u,
     +        kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,
     +        kimax,kimax2,maxiter,imix,igrd,ndvgrd,idsprs,isprsv,
     +        idsprs0,idsprsl,idsprsi,idsprsv,kk,layers,nnne
      INTEGER lepr,nwd,irecl,irecl0,irecl1,lmd,lmplmd,iwd,ikpt,
     +        gw,gw_neigd,ir_fac,int_len,nkpt_l
      REAL alpha,chng,delgau,e1s,e2s,ef,gmaxxc,seigc,
     +     seigscv,sprsv,start_it,time1,time2,timeq1,timeq2,tkb,
     +     ts,te_vcoul,te_veff,te_exc,spinf,e2_dos,e1_dos,sig_dos
      REAL area,delz,dvac,gmax,omtil,vol,volint,gmax_inp,
     +     z1,zksum,theta,phi,sigma,zsigma,e_u_c,scale
      LOGICAL film,invs,invs2,pot8,symor,total,vchk,zrfs,l_wann,stop80 
      LOGICAL l_f,ctail,disp,dos,eonly,form66,frcor,gauss,
     +        iplot,lflip,plpot,reap,score,secvar,slice,strho,swsp,
     +        tria,lwb,cdinf,l_endit,l_f2u,l_u2f,l_bmt,l_zref,form76,
     +        vacdos,integ,pallst,l_soc,l_opti,dos_tmp
      CHARACTER(len=3) latnam
      CHARACTER(len=4) namgrp
      CHARACTER(len=9) ivers
      CHARACTER(len=12) fname(3)
c
c+dw new arguments for STM
      INTEGER nstars,nstm
      REAL    tworkf,locx(2),locy(2)
      LOGICAL starcoeff

c---> pk non-collinear
c---> new variables for non-collinear calculations
      REAL qss(3),mix_b
      LOGICAL l_noco,l_ss,l_mperp,l_constr
C     UNIT NUMBERS
      INTEGER ndomfile,npotmatfile,noinpfile,nrhomfile
c---> pk non-collinear

c--- J<
      INTEGER             :: nqpt,phn,phnd,nsh,mtypes,itype
      REAL                :: thetaJ,qn
      REAL, PARAMETER     :: tol = 1.e-8 
      INTEGER             :: qcount,nmagn,imt,i_J,j_J
      LOGICAL             :: l_J,l_disp,l_wr,l_jenerg
      LOGICAL, ALLOCATABLE:: l_magn(:)
      REAL, ALLOCATABLE   :: M(:)
      INTEGER, ALLOCATABLE:: magtype(:),nmagtype(:)
c--- J>

C     ..
C     .. Local Arrays ..
      REAL amat(3,3),bmat(3,3),bbmat(3,3),sig_b(2)
      LOGICAL eig66(2),eigvar(3), evonly(2)
      LOGICAL sso_opt(2)
      CHARACTER*2 namat(0:103)
      CHARACTER*8 name(10)
c
      COMPLEX, ALLOCATABLE :: clnu(:,:,:)
      COMPLEX, ALLOCATABLE :: d_wgn(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ustep(:)
      INTEGER, ALLOCATABLE :: ig(:,:,:)
      INTEGER, ALLOCATABLE :: ig2(:)
      INTEGER, ALLOCATABLE :: igz(:)
      INTEGER, ALLOCATABLE :: jri(:)
      INTEGER, ALLOCATABLE :: kv2(:,:)
      INTEGER, ALLOCATABLE :: kv3(:,:)
      INTEGER, ALLOCATABLE :: llh(:,:)
      INTEGER, ALLOCATABLE :: lmax(:)
      INTEGER, ALLOCATABLE :: mlh(:,:,:)
      INTEGER, ALLOCATABLE :: mrot(:,:,:)
      INTEGER, ALLOCATABLE :: ncv(:)
      INTEGER, ALLOCATABLE :: neq(:)
      INTEGER, ALLOCATABLE :: ngopr(:)
      INTEGER, ALLOCATABLE :: nlh(:)
      INTEGER, ALLOCATABLE :: nmem(:,:)
      INTEGER, ALLOCATABLE :: nstr(:)
      INTEGER, ALLOCATABLE :: nstr2(:)
      INTEGER, ALLOCATABLE :: ntypsy(:)
      INTEGER, ALLOCATABLE :: igfft(:,:)
      INTEGER, ALLOCATABLE :: igfft2(:,:)
      INTEGER, ALLOCATABLE :: nflip(:)
      INTEGER, ALLOCATABLE :: kq1_fft(:)
      INTEGER, ALLOCATABLE :: kq2_fft(:)
      INTEGER, ALLOCATABLE :: kq3_fft(:)
      INTEGER, ALLOCATABLE :: nq3_fft(:)
      INTEGER, ALLOCATABLE :: kmxq_fft(:)
      INTEGER, ALLOCATABLE :: ncst(:)
      INTEGER, ALLOCATABLE :: izlay(:,:)
      INTEGER, ALLOCATABLE :: invarop(:,:)
      INTEGER, ALLOCATABLE :: invarind(:)
      INTEGER, ALLOCATABLE :: multab(:,:)
      INTEGER, ALLOCATABLE :: invtab(:)
      INTEGER, ALLOCATABLE :: invsat(:)
      INTEGER, ALLOCATABLE :: invsatnr(:)
      INTEGER, ALLOCATABLE :: lnonsph(:)
      INTEGER, ALLOCATABLE :: nkpt(:)
      INTEGER, ALLOCATABLE :: nlo(:)
      INTEGER, ALLOCATABLE :: llo(:,:)
      INTEGER, ALLOCATABLE :: nlol(:,:)
      INTEGER, ALLOCATABLE :: lo1l(:,:)
      INTEGER, ALLOCATABLE :: skiplo(:,:)
      INTEGER, ALLOCATABLE :: lapw_l(:)
      INTEGER, ALLOCATABLE :: ulo_der(:,:)
      REAL, ALLOCATABLE :: dx(:)
      REAL, ALLOCATABLE :: pos(:,:)
      REAL, ALLOCATABLE :: rgphs(:,:,:)
      REAL, ALLOCATABLE :: rmsh(:,:)
      REAL, ALLOCATABLE :: rmt(:)
      REAL, ALLOCATABLE :: sk2(:)
      REAL, ALLOCATABLE :: phi2(:)
      REAL, ALLOCATABLE :: sk3(:)
      REAL, ALLOCATABLE :: tau(:,:)
      REAL, ALLOCATABLE :: taual(:,:)
      REAL, ALLOCATABLE :: volmts(:)
      REAL, ALLOCATABLE :: zatom(:)
      REAL, ALLOCATABLE :: el0(:,:,:,:)
      REAL, ALLOCATABLE :: evac0(:,:,:)
      REAL, ALLOCATABLE :: ello0(:,:,:)
      REAL, ALLOCATABLE :: force(:,:,:)
      REAL, ALLOCATABLE :: force_old(:,:)
      REAL, ALLOCATABLE :: bk(:,:,:)
      REAL, ALLOCATABLE :: wtkpt(:,:)
      REAL, ALLOCATABLE :: ellow(:)
      REAL, ALLOCATABLE :: elup(:)
      REAL, ALLOCATABLE :: rkmax(:)
      REAL, ALLOCATABLE :: zelec(:)
      REAL, ALLOCATABLE :: pgfft(:)
      REAL, ALLOCATABLE :: pgfft2(:)
      REAL, ALLOCATABLE :: ufft(:)
      REAL, ALLOCATABLE :: bmu(:)
      REAL, ALLOCATABLE :: vr0(:)
      REAL, ALLOCATABLE :: w_iks(:,:,:)
      REAL, ALLOCATABLE :: pgft2x(:)
      REAL, ALLOCATABLE :: pgft2y(:)
      REAL, ALLOCATABLE :: pgft2xx(:)
      REAL, ALLOCATABLE :: pgft2yy(:)
      REAL, ALLOCATABLE :: pgft2xy(:)
      REAL, ALLOCATABLE :: enmix(:,:) 
      REAL, ALLOCATABLE :: qj(:,:)
      REAL, ALLOCATABLE :: eig_l(:,:)
      LOGICAL, ALLOCATABLE :: lchange(:,:,:,:),l_dulo(:,:)
      LOGICAL, ALLOCATABLE :: lchg_v(:,:,:),l_geo(:)
      LOGICAL, ALLOCATABLE :: llochg(:,:,:),l_relax(:)
      LOGICAL, ALLOCATABLE :: soc_opt(:) 
c---> pk non-collinear
c---> new arrays for non-collinear calculations
      INTEGER, ALLOCATABLE :: nq2_fft(:)
      INTEGER, ALLOCATABLE :: kmxq2_fft(:)
      REAL, ALLOCATABLE :: alph(:),beta(:),b_con(:,:),alph1(:)
c---> pk non-collinear
c+lda+u
      TYPE (t_utype), ALLOCATABLE :: lda_u(:)
c-lda+u
c+odim
      INTEGER, TARGET, ALLOCATABLE :: ig1(:,:)
      INTEGER, TARGET, ALLOCATABLE :: kv1(:,:)
      INTEGER, TARGET, ALLOCATABLE :: nstr1(:)
      INTEGER, TARGET, ALLOCATABLE :: ngopr1(:)
      REAL,    TARGET, ALLOCATABLE :: mrot1(:,:,:)
      REAL,    TARGET, ALLOCATABLE :: tau1(:,:)
      INTEGER, TARGET, ALLOCATABLE :: invtab1(:)
      INTEGER, TARGET, ALLOCATABLE :: multab1(:,:)
      INTEGER, TARGET, ALLOCATABLE :: igfft1(:,:)
      REAL,    TARGET, ALLOCATABLE :: pgfft1(:)
      REAL,    TARGET, ALLOCATABLE :: pgft1x(:)
      REAL,    TARGET, ALLOCATABLE :: pgft1y(:)
      REAL,    TARGET, ALLOCATABLE :: pgft1xx(:)
      REAL,    TARGET, ALLOCATABLE :: pgft1yy(:)
      REAL,    TARGET, ALLOCATABLE :: pgft1xy(:)

      TYPE (od_dim) :: odd
      TYPE (od_inp) :: odi
      TYPE (od_sym) :: ods
      TYPE (od_lda) :: odl
      TYPE (od_gga) :: odg
c-odim
C     ..
C     .. DATA statements ..
      DATA namat/'va',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     +     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     +     ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     +     'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     +     'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce',
     +     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +     'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     +     'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu',
     +     'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw'/
C     ..
c+t3e
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER isize,irank,taskID
      INTEGER ierr(3)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,isize,ierr)
#else
      INTEGER isize,irank
      irank=0 ; isize=1
#endif
!
!open hdf-library
!
#ifdef DUMP_DATA_HDF
! open hdf5, probably conflicts with CPP_HDF, which is not used
      call binmat_start()
#endif
#ifdef CPP_HDF
#ifdef CPP_HPM
      call f_hpminit( taskID, "fleur" ) 
      call f_hpmstart( 1, "total it"  ) 
#endif
      CALL hdf_init()
#endif
!
! names of direct access-files
!
#if ( defined(CPP_MPI) && !defined(CPP_HDF) )
      aoff = iachar('1')-1
      d1  = mod(irank,10)
      IF (irank < 100) THEN
        d10 = int( (irank + 0.5)/10 )
        fname(1) = 'eig'//achar(d10+aoff)//achar(d1+aoff)
      ELSE
        d10 = mod((irank-d1)/10,10)
        d100 = (irank-10*d10-d1)/100
        IF ( d100.GE.10 ) d100 = d100 + iachar('7')
        fname(1) =
     +  'eig'//achar(d100+aoff)//achar(d10+aoff)//achar(d1+aoff)
      ENDIF
#else
      fname(1) = 'eig'
#endif
      fname(2) = 'tmat'
      fname(3) = 'tmas'

      ivers='fleur.26'
      seigscv = 0.0
      te_vcoul = 0.0
      te_veff = 0.0
      te_exc = 0.0
      e_u_c = 0.0
      ts = 0.0
      gw       = -1
      gw_neigd =  0
c-t3e
      IF (irank.EQ.0) THEN
        OPEN (2,file='time.info',form='formatted',status='unknown')
        OPEN (6,file='out',form='formatted',status='unknown')
      ENDIF

      CALL dimens(
     >            irank,ivers,gw,gw_neigd,
     <            nop,k1d,k2d,k3d,n3d,n2d,
     <            kq1d,kq2d,kq3d,kxc1d,kxc2d,kxc3d,
     <            ntypd,natd,jmtd,ntypsd,nlhd,memd,
     <            lmaxd,jspd,nvacd,nmzd,nmzxyd,nspd,
     <            nvd,nv2d,nwdd,nkptd,
     <            neigd,msh,nstd,ncvd,layerd,nlod,llod,odd,
     <            itmax)

      IF (irank.EQ.0) THEN
#ifdef CPP_HTML
        OPEN (16,file='inf.html',form='formatted',status='unknown')
        CALL html_out('top',0,itmax,ivers)
#else
        OPEN (16,file='inf',form='formatted',status='unknown')
#endif
      ENDIF

      nn2d= (2*k1d+1)* (2*k2d+1)
      nn3d= (2*k1d+1)* (2*k2d+1)* (2*k3d+1)
c-odim
      IF (odd%d1) THEN
         odd%k3 = k3d
         odd%nn2d = (2*(odd%k3) + 1)*(2*(odd%M) + 1)
      ELSE
         odd%k3 = 0 ; odd%M =0 ; odd%nn2d = 1
      ENDIF
c-odim
      ALLOCATE ( clnu(memd,0:nlhd,ntypsd),ustep(n3d) )
      ALLOCATE ( ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d),igz(n3d) )
      ALLOCATE ( jri(ntypd),kv2(2,n2d),kv3(3,n3d),llh(0:nlhd,ntypsd) )
      ALLOCATE ( lmax(ntypd),mlh(memd,0:nlhd,ntypsd),mrot(3,3,nop) )
      ALLOCATE ( ncv(ntypd),neq(ntypd),ngopr(natd) )
      ALLOCATE ( nlh(ntypsd),nmem(0:nlhd,ntypsd) )
      ALLOCATE ( nstr2(n2d),ntypsy(natd),nstr(n3d) )
      ALLOCATE ( igfft(0:nn3d-1,2),igfft2(0:nn2d-1,2),nflip(ntypd) )
      ALLOCATE ( kq1_fft(nwdd),kq2_fft(nwdd),kq3_fft(nwdd) )
      ALLOCATE ( nq2_fft(nwdd),nq3_fft(nwdd) )
      ALLOCATE ( kmxq2_fft(nwdd),kmxq_fft(nwdd),ncst(ntypd) )
      ALLOCATE ( izlay(layerd,2) )
      ALLOCATE ( invarop(natd,nop),invarind(natd) )
      ALLOCATE ( multab(nop,nop),invtab(nop) )
      ALLOCATE ( invsat(natd),invsatnr(natd) )
      ALLOCATE ( lnonsph(ntypd),nkpt(nwdd) )
      ALLOCATE ( dx(ntypd),pos(3,natd),tau(3,nop) )
      ALLOCATE ( rmsh(jmtd,ntypd),rmt(ntypd),sk2(n2d),sk3(n3d) )
      ALLOCATE ( phi2(n2d) )
      ALLOCATE ( taual(3,natd),volmts(ntypd),zatom(ntypd) )
      ALLOCATE ( el0(0:lmaxd,ntypd,jspd,nwdd) )
      ALLOCATE ( evac0(2,jspd,nwdd),rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)  )
      ALLOCATE ( force(3,ntypd,jspd) )
      ALLOCATE ( force_old(3,ntypd) )
      ALLOCATE ( bk(3,nkptd,nwdd),wtkpt(nkptd,nwdd) )
      ALLOCATE ( ellow(nwdd),elup(nwdd),rkmax(nwdd) )
      ALLOCATE ( zelec(nwdd) )
      ALLOCATE ( pgfft(0:nn3d-1),pgfft2(0:nn2d-1) )
      ALLOCATE ( ufft(0:27*k1d*k2d*k3d-1) )
      ALLOCATE ( bmu(ntypd),vr0(ntypd) )
      ALLOCATE ( lchange(0:lmaxd,ntypd,jspd,nwdd) )
      ALLOCATE ( lchg_v(2,jspd,nwdd),l_geo(ntypd) )
      ALLOCATE ( nlo(ntypd),llo(nlod,ntypd),skiplo(ntypd,jspd) )
      ALLOCATE ( ello0(nlod,ntypd,jspd),llochg(nlod,ntypd,jspd) )
      ALLOCATE ( lo1l(0:llod,ntypd),nlol(0:llod,ntypd),lapw_l(ntypd) )
      ALLOCATE ( alph(ntypd),beta(ntypd),l_relax(ntypd) )
      ALLOCATE ( alph1(ntypd),l_magn(ntypd),M(ntypd) )
      ALLOCATE ( magtype(ntypd),nmagtype(ntypd) )
      ALLOCATE ( b_con(2,ntypd),lda_u(ntypd),l_dulo(nlod,ntypd) )
      ALLOCATE ( enmix(jspd,nwdd),d_wgn(-3:3,-3:3,3,nop) )
      ALLOCATE ( ulo_der(nlod,ntypd) )
      ALLOCATE ( soc_opt(ntypd+2) )
c+odim
      ALLOCATE ( ig1(-odd%k3:odd%k3,-odd%M:odd%M) )
      ALLOCATE ( kv1(2,odd%n2d),nstr1(odd%n2d) )
      ALLOCATE ( ngopr1(natd),mrot1(3,3,odd%nop),tau1(3,odd%nop) )
      ALLOCATE ( invtab1(odd%nop),multab1(odd%nop,odd%nop) )
      ALLOCATE ( igfft1(0:odd%nn2d-1,2),pgfft1(0:odd%nn2d-1) )
      sk2(:) = 0.0 ; phi2(:) = 0.0
c-odim

c---> pk non-collinear
      noinpfile   = 24
      ndomfile    = 22
      nrhomfile   = 26
      npotmatfile = 25
c---> pk non-collinear

      DO n = 1,ntypd
        vr0(n) = 0.0
        M(n)= 0.0
        l_magn(n)=.false.
        DO jspin=1,jspd
          DO i=1,3
            force(i,n,jspin) = 0.0
          ENDDO
        ENDDO
      ENDDO

      time1 = 0.0
      time2 = 0.0
      start_it = 0.0
      CALL cpu_time(time1)
c--- J< 
      timeq1 = 0.0
      timeq2 = 0.0
      l_wr=.true.
      nqptd=1
      nmagn=1
      mtypes=1
      phnd=1
c--- J>
c+t3e
      IF (irank.EQ.0) THEN
c-t3e
      CALL inped( 
     >           ntypd,natd,nwdd,layerd,nmzxyd,nmzd,jmtd,lmaxd,
     >           nlod,llod,ivers,namat,noinpfile,
     <           strho,dos,isec1,ndir,secvar,icorr,krla,namgrp,
     <           cdinf,ncst,pot8,vchk,latnam,name,
     <           igrd,ndvgrd,idsprs,isprsv,idsprs0,idsprsl,
     <           idsprsi,idsprsv,chng,sprsv,lwb,gmaxxc,
     <           gauss,tria,delgau,tkb,ellow,elup,nwd,zelec,
     <           lpr,form66,l_f,eonly,eig66,lnonsph,rkmax,
     <           frcor,slice,ctail,disp,kcrel,lepr,
     <           itmax,maxiter,imix,alpha,spinf,l_u2f,l_f2u,l_bmt,
     <           film,invs,invs2,zrfs,omtil,vol,volint,volmts,
     <           area,z1,delz,lmax,l_soc,soc_opt,theta,phi,l_geo,
     <           jspins,nvac,nmz,nmzxy,dvac,scale,ntype,neq,jri,gmax,
     <           swsp,lflip,bmu,nflip,vacdos,layers,integ,
     <           izlay,iplot,score,plpot,kk,e1s,e2s,nnne,pallst,
     <           amat,bmat,dx,pos,taual,rmt,rmsh,zatom,lapw_l,ulo_der,
     <           nlo,llo,nlol,lo1l,l_dulo,tworkf,nstars,nstm,starcoeff,
     <           locx,locy,l_noco,l_ss,l_mperp,l_constr,mix_b,alph,beta,
     <           b_con,l_relax,qss,e1_dos,e2_dos,sig_dos,lda_u,n_u,
     <           gw,gw_neigd,odd,
     <           l_J,thetaJ,l_magn,nmagn,M,mtypes,
     <           magtype,nmagtype,nsh,l_disp,sso_opt)
c
      IF (igrd.NE.0) THEN
         ALLOCATE (pgft2x(0:nn2d-1),pgft2xx(0:nn2d-1),pgft2xy(0:nn2d-1),
     +             pgft2y(0:nn2d-1),pgft2yy(0:nn2d-1))
c-odim
         ALLOCATE (pgft1x(0:odd%nn2d-1),pgft1xx(0:odd%nn2d-1),
     +             pgft1xy(0:odd%nn2d-1),
     +             pgft1y(0:odd%nn2d-1),pgft1yy(0:odd%nn2d-1))
      ELSE
         ALLOCATE (pgft2x(0:1),pgft2xx(0:1),pgft2xy(0:1),
     +             pgft2y(0:1),pgft2yy(0:1))
c-odim
         ALLOCATE (pgft1x(0:1),pgft1xx(0:1),pgft1xy(0:1),
     +             pgft1y(0:1),pgft1yy(0:1))
      ENDIF
      odd%nq2 = odd%n2d
c-odim
c+t3e
      IF (.not.(strho.OR.l_f2u.OR.l_u2f.OR.iplot)) THEN
        IF (l_noco) THEN
          INQUIRE (file='rhomat_inp',exist=strho) ! if no density (rhomat_inp) exists,
        ELSE
          INQUIRE (file='cdn1',exist=strho)       ! if no density (cdn1) exists,
        ENDIF
        strho = .not.strho                ! create a starting density
      ENDIF
      l_opti = .false.
      IF ((iplot).OR.(strho).OR.(swsp).OR.
     +    (lflip).OR.(l_f2u).OR.(l_u2f).OR.(l_bmt)) l_opti = .true.
        
      form76 = .false.
      IF (l_soc.AND.form66) THEN
        IF (.NOT.eonly) STOP 'SOC: form66=T only with eonly=T !'
        form66 = .false.
        form76 = .true.
      ENDIF
c
      CALL setup(
     >     odd,lmaxd,nkptd,jspd,ncvd,natd,ntypd,nlhd,memd,nstd,
     >     ntypsd,nwdd,neigd,nop,nn3d,nn2d,n2d,n3d,k3d,k2d,k1d,
     >     kxc3d,kxc2d,kxc1d,kq3d,kq2d,kq1d,taual,ncst,zelec,phi,
     >     nmzxy,nmz,jmtd,nmzxyd,nmzd,invs2,delz,n_u,l_soc,theta,
     >     invs,film,zrfs,l_f,area,dvac,scale,z1,vol,omtil,l_opti,
     >     igrd,jspins,l_ss,ntype,nwd,nvac,rkmax,gmax,gmaxxc,
     >     neq,lmax,amat,bmat,iplot,latnam,namgrp,volmts,rmt,gw,
     <     bbmat,llh,mlh,nmem,clnu,nlh,nsymt,ntypsy,
     <     invarind,invarop,invtab,multab,invsatnr,invsat,
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
#ifdef DUMP_DATA_OLD
! output amat, bmat
      write(misc_name, '(A, ".bin")') 'real_basis'
      call write_mat(amat, misc_name)
      write(misc_name, '(A, ".bin")') 'rec_basis'
      call write_mat(bmat, misc_name)
! atom positions
      write(misc_name, '(A, ".bin")') 'positions'
      call write_mat(taual, misc_name)
! TODO: scalar file, need: lmax, num_atoms, max_gpoints (nvd)
! and for the radmesh: energies(0..lmax), jri,dx,rmt for each atom
      open(321, file='scalars.dat', status='replace')
      write(321, *) lmax, natd, nvd
! mesh parameters are only stored per atom type. that's not in the code (yet),
! so just write them out multiple times
      do i = 1, natd
        write(321, *) jri(1), dx(1), rmt(1)
      end do
! energies are not available yet
!      write(321, *) el
      close(321)
#endif
#ifdef DUMP_DATA_HDF
! HDF output:
! * Basis matrices
      call write_hdf_mat('real_basis', amat)
      call write_hdf_mat('rec_basis', bmat)
! * Atoms: coordinates, type and rotation matrix (identity if atom is represantative)
      call write_hdf_atoms(natd, taual, mrot(:, :, invtab(ngopr(:))),
     +                      ntype, neq)
! * l-cutoffs
      call write_hdf_mat('lmax', lmax,
     +   "l-cutoff, one entry per atom type")
      call write_hdf_mat('lmax_nonspherical', lnonsph,
     +   "l-cutoff for the 'lower' part @ hssphn, " //
     +    "one entry per atom type")
! * mesh parameters
      call write_hdf_mesh(ntype, jri, dx, rmt)
#endif
c
      ng3=nq3 ; ng2=nq2 
      CALL cpu_time(time2)
      WRITE (2,FMT='(/,t20,a,/)') 'FLAPW7 time.info-file'
      CALL outtime('preparation: generation of stars,lattice harmonics,
     +etc:',time2-time1)
c+t3e
      ENDIF ! irank.eq.0
c
c-odim
      odd%nq2 = odd%n2d
      odd%kimax2 = odd%nq2 - 1
      odd%nat = natd

      odi%d1 = odd%d1 ; odi%mb = odd%mb ; odi%M = odd%M
      odi%k3 = odd%k3 ; odi%chi = odd%chi ; odi%rot = odd%rot
      odi%invs = odd%invs ; odi%zrfs = odd%zrfs
      odi%n2d = odd%n2d ; odi%nq2 = odd%nq2 ; odi%nn2d = odd%nn2d
      odi%kimax2 = odd%kimax2 ; odi%m_cyl = odd%m_cyl
      odi%ig => ig1 ; odi%kv => kv1 ; odi%nst2 => nstr1

      ods%nop = odd%nop ; ods%nat = odd%nat
      ods%mrot => mrot1 ; ods%tau => tau1 ; ods%ngopr => ngopr1
      ods%invtab => invtab1 ; ods%multab => multab1

      odl%nn2d = odd%nn2d
      odl%igf => igfft1 ; odl%pgf => pgfft1

      odg%nn2d = odd%nn2d
      odg%pgfx => pgft1x ; odg%pgfy => pgft1y
      odg%pgfxx => pgft1xx ; odg%pgfyy => pgft1yy ; odg%pgfxy => pgft1xy
c+odim
c
#ifdef CPP_MPI
      CALL MPI_BCAST(l_opti,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_noco,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_soc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(strho ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jspins,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      IF (l_opti) THEN
         IF (iplot .AND. (irank==0) ) THEN
            IF (l_noco) THEN
            CALL pldngen(
     >                  jspd,nop,k1d,k2d,k3d,n2d,n3d,
     >                  ntypd,ntypsd,nlhd,jmtd,natd,nmzd,nmzxyd,
     >                  nrhomfile,neq,invs,invs2,rmt,rmsh,zatom,namat,
     >                  dx,delz,z1,sk3,tau,bmat,ig,ig2,symor,vol,taual,
     >                  jspins,film,nvac,volmts,volint,area,
     >                  ntype,jri,ntypsy,nlh,alph,beta,
     >                  mrot,ng2,ng3,kv3,igfft2,invtab,
     >                  igfft,pgfft2,pgfft,nstr2,nstr,
     >                  kimax,kimax2,nmz,nmzxy,sigma,odi,slice)
           ENDIF
         ENDIF
         CALL optional(
     >    irank,isize,
     >    ntypd,ntypsd,nmzxyd,nmzd,nlhd,jmtd,msh,memd,
     >    jspd,n3d,n2d,ntype,jri,jspins,lmaxd,nop,natd,
     >    k1d,k2d,k3d,nspd,nstd,mx3,omtil,zrfs,ig,taual,
     >    iplot,strho,swsp,lflip,l_f2u,l_u2f,l_bmt,symor,sk3,
     >    plpot,idsprs,ndvgrd,igrd,krla,icorr,lepr,
     >    name,namat,nq3,nq2,nmzxy,nmz,delz,z1,volint,
     >    invs2,invs,film,nflip,bmat,amat,pos,invtab,
     >    neq,nstr2,nstr,lmax,mlh,llh,nmem,clnu,ngopr,
     >    ntypsy,nlh,slice,bmu,sprsv,isprsv,mrot,tau,
     >    nvac,zatom,dx,rmt,rmsh,kv2,kv3,nop2,
     >    ig2,vol,volmts,area,score,l_noco,nwd,nlod,nlo,
     >    llo,l_dulo,nwdd,ellow,elup,sigma,lda_u(:)%l,n_u,
     >    sk2,phi2,odi,ods)
      ENDIF
c
c--> determine recl for direct access files (28,29,38 and 66)
c
      nbasfcn = nvd + natd*nlod*(2*llod+1)
      irecl = ntypd * ( lmaxd+1 + nlod) + 2
c---> in non-collinear calculations, the energy parameters have to be 
c---> saved for both spins.
      IF (l_noco) irecl = 2*irecl
c
c IF no-inversion, the eigenvectors are COMPLEX (2*nvd)
c IF spin-orbit coupling, they are COMPLEX and twice as large
c the 5 is ONLY needed for spin-spiral, because nv depends on jsp
      IF (l_soc.OR.l_noco) THEN
        irecl = 8 * (irecl + 5 + 2*neigd*(2*nbasfcn + 1))
!       IF(l_J) THEN
!        irecl = 8 * (irecl + 5 + 2*neigd)
!       ENDIF
      ELSE
#ifdef CPP_INVERSION
        irecl = 8 * (irecl + 4 + neigd*(nbasfcn + 1))
!        IF (n_u.GT.0) irecl = irecl + 8 * neigd*nbasfcn
#else
        irecl = 8 * (irecl + 4 + neigd*(2*nbasfcn + 1))
#endif
      ENDIF
!
! --> add some information for the K-vectors fo the LO's
!
      irecl = irecl + 4 * natd*nlod*(2*llod+1)
c
      lmd     = lmaxd* (lmaxd+2)
      lmplmd  = (lmd* (lmd+3))/2
      irecl0 = 8* (8* (lmplmd+1)+ (lmd+1)**2)
      irecl0 = irecl0 + 16*(2*llod+1)*( nlod*2*(lmd+1) + 
     +                     (2*llod+1)*nlod*(nlod+1)/2 )
#ifdef CPP_MPI
         CALL MPI_BCAST(n_u,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      IF (n_u > 0) irecl0 = irecl0 + 16*49
c+gu
      l_zref=.false.
      IF ((irank.EQ.0).AND.(.NOT.l_opti)) THEN
        zksum=0.0
        DO iwd=1,nwd
          DO ikpt=1,nkpt(iwd)
            zksum=zksum+ABS(bk(3,ikpt,iwd))
          ENDDO
        ENDDO
        IF (zrfs.AND.(zksum.LT.1e-9)) l_zref=.true.
      ENDIF
!      IF (isize > 1)  l_zref = .false.
      IF (irank.EQ.0) THEN
!         DO n = 1,ntype
!            IF (nlo(n)>=1) l_zref = .false.
!         ENDDO
         IF (l_noco) l_zref = .FALSE.

         ir_fac = size(transfer(zksum,(/1/)))
         IF ( ir_fac == 1 ) THEN
             int_len = 8                 ! ir_fac = 1 if integer*8
         ELSEIF ( ir_fac == 2 ) THEN
             int_len = 4                 ! ir_fac = 2 if integer*4
         ELSE
             STOP 'fleur: strange length of integers'
         ENDIF
         irecl = irecl + int_len *( 3 + 3*nvd )
         IF (l_ss) irecl = irecl + int_len *( 3 + 2*3*nvd )
         IF (l_J) THEN
           irecl = 8*(neigd+3) + 4*int_len
         ENDIF
         IF (l_zref) THEN
            irecl1 = (5*nvd+2)*int_len + (nvd+3)*8
         ELSE
            irecl1 = (7*nvd+8)*int_len - (2*nvd+1)*int_len
         ENDIF

         IF (l_ss) irecl1 = 2*irecl1
c
c--- J< 
        l_jenerg = .false.
        IF (l_J) THEN
         OPEN (113,file='qpts',status='old')
         INQUIRE(file='jenerg',exist=l_jenerg)
         OPEN (114,file='jenerg',status='unknown')
         READ (113,*) nqptd
        ENDIF
      ENDIF!(irank.EQ.0)

#ifdef CPP_MPI
         CALL MPI_BCAST(nqptd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(l_jenerg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif
        ALLOCATE ( qj(3,nqptd) )
        nqpt=nqptd
c+t3e
      IF (irank.EQ.0) THEN
c-t3e
       IF (l_J) THEN
        IF(l_disp)THEN
        WRITE(6,*) 'q points for the magnon spectrum calculation:'
        ELSE
        WRITE(6,*) 'q points for the J-constants calculation:'
        ENDIF
        WRITE(6,*)'      q1       ','      q2       ','      q3'
        DO qcount=1,nqpt
          READ (113,3333) qj(1,qcount),qj(2,qcount),qj(3,qcount)
          WRITE(6,3333) qj(1,qcount),qj(2,qcount),qj(3,qcount)
        ENDDO
 3333   FORMAT(3(f14.10,1x))
       ENDIF !(l_J)
c+t3e
      ENDIF
c-t3e
c--- J>


#ifdef CPP_MPI
      CALL mpi_bc_all(
     >           n3d,memd,nlhd,ntypsd,ntypd,nwdd,natd,n2d,
     >           nop,k1d,k2d,k3d,nkptd,nqptd,jspd,lmaxd,jmtd,nlod,
     >           irecl,irecl0,irecl1,l_zref,nop2,llod,
     X           ntype,lpr,form66,eig66,eonly,secvar,isec1,
     X           l_f,cdinf,dos,vacdos,integ,pallst,nn3d,
     X           iplot,strho,swsp,lflip,l_f2u,l_u2f,l_bmt,slice,
     X           zrfs,film,invs,invs2,jspins,nvac,neq,lmax,
     X           bk,wtkpt,rkmax,ellow,elup,nwd,nkpt,nsymt,
     X           jri,rmt,dx,ntypsy,nlh,nq2,nq3,nmz,nmzxy,
     X           ustep,bbmat,bmat,invsat,ncv,pot8,invsatnr,zatom,
     X           nstr,nstr2,kv2,kv3,sk3,igfft,pgfft,kmxq_fft,nq3_fft,
     X           kq1_fft,kq2_fft,kq3_fft,zelec,izlay,ello0,
     X           ig2,ig,rgphs,area,delz,z1,omtil,volint,vol,volmts,
     X           e1s,e2s,kk,layers,nnne,ndir,mx1,mx2,mx3,layerd,
     X           lnonsph,el0,evac0,lepr,irank,itmax,alpha,nmem,mlh,llh,
     X           rmsh,clnu,ngopr,mrot,tau,taual,invarop,invarind,multab,
     X           invtab,amat,ncst,d_wgn,l_soc,soc_opt,theta,phi,
     X           nlo,llo,nlol,lo1l,l_dulo,ulo_der,ctail,skiplo,tworkf,
     X           nstars,nstm,starcoeff,locx,locy,l_noco,l_ss,l_mperp,
     X           l_constr,alph,beta,b_con,l_relax,qss,lda_u,igrd,lapw_l,
     X           mix_b,kmxq2_fft,nq2_fft,l_geo,neigd,sigma,zsigma,sig_b,
     X           pos,sk2,phi2,odd,odd%M,odd%nn2d,odd%n2d,odd%nop,igfft1,
     X           pgfft1,ig1,kv1,nstr1,ngopr1,tau1,invtab1,multab1,mrot1,
     X           l_J,l_disp,sso_opt,qj,nmagn,mtypes,gw,gw_neigd,
     X           l_magn,thetaJ,magtype,nmagtype,
     <           TYP_REAL,TYP_COMPLEX)

      IF (irank >  0) THEN 
        IF (igrd.NE.0) THEN
         ALLOCATE (pgft2x(0:nn2d-1),pgft2xx(0:nn2d-1),pgft2xy(0:nn2d-1),
     +             pgft2y(0:nn2d-1),pgft2yy(0:nn2d-1))
        ELSE
         ALLOCATE (pgft2x(0:1),pgft2xx(0:1),pgft2xy(0:1),
     +             pgft2y(0:1),pgft2yy(0:1))
        ENDIF
      ENDIF

#else
      TYP_REAL= 0
#endif 
c
      IF (iplot) STOP 'density plot o.k.'
      IF (strho) STOP 'starting density generated'
      IF (swsp)  STOP 'spin polarised density generated'
      IF (lflip) STOP 'magnetic moments flipped'
      IF (l_f2u) STOP 'conversion to unformatted'
      IF (l_u2f) STOP 'conversion to formatted'
      IF (l_bmt) STOP '"cdnbmt" written'
      qn = 0.0
c-t3e
c-odim
      odd%nq2 = odd%n2d
      odd%kimax2 = odd%nq2 - 1
      odd%nat = natd

      odi%d1 = odd%d1 ; odi%mb = odd%mb ; odi%M = odd%M
      odi%k3 = odd%k3 ; odi%chi = odd%chi ; odi%rot = odd%rot
      odi%invs = odd%invs ; odi%zrfs = odd%zrfs
      odi%n2d = odd%n2d ; odi%nq2 = odd%nq2 ; odi%nn2d = odd%nn2d
      odi%kimax2 = odd%kimax2 ; odi%m_cyl = odd%m_cyl
      odi%ig => ig1 ; odi%kv => kv1 ; odi%nst2 => nstr1

      ods%nop = odd%nop ; ods%nat = odd%nat
      ods%mrot => mrot1 ; ods%tau => tau1 ; ods%ngopr => ngopr1
      ods%invtab => invtab1 ; ods%multab => multab1

      odl%nn2d = odd%nn2d
      odl%igf => igfft1 ; odl%pgf => pgfft1

      odg%nn2d = odd%nn2d
      odg%pgfx => pgft1x ; odg%pgfy => pgft1y
      odg%pgfxx => pgft1xx ; odg%pgfyy => pgft1yy ; odg%pgfxy => pgft1xy
c+odim
      IF (l_noco) nbasfcn = 2*nbasfcn
c
c--- J<
      IF (l_J) THEN
        itmax = 1
        phnd=2
        nkpt_l = CEILING(REAL(nkptd)/isize)
        ALLOCATE ( eig_l(neigd+5,nkpt_l) )
      ELSE
        nkpt_l = 1
        ALLOCATE ( eig_l(1,nkpt_l) )
      ENDIF
c--- J>

      DO 80 it = 1,itmax
c+t3e
         IF (alpha.LT.10.0) THEN

         IF (it.GT.1) THEN
            eig66(1)= .false. 
         END IF
         IF (irank.EQ.0) THEN
c-t3e
           CALL cpu_time(start_it)
           WRITE (2,FMT='(/,a,i3,/)') 'ITERATION',it
c
           IF (it.GT.1) THEN
              pot8 = .false.
              alpha = alpha - NINT(alpha)
           END IF
c
           WRITE (6,FMT=8100) it
           WRITE (16,FMT=8100) it
 8100      FORMAT (/,10x,'   it=    ',i5)
c
           IF (.not.pot8) THEN
c
c      ----> potential generator
c
c---> pk non-collinear
c--->        reload the density matrix from file rhomat_in
c--->        calculate spin-up and -down density for USE in the
c--->        potential generator and store the direction of 
c--->        magnetization on file dirofmag
             IF (l_noco) THEN
               CALL cpu_time(time1)
               CALL rhodirgen(
     >                     jspd,nop,k1d,k2d,k3d,n2d,n3d,ntypd,ntypsd,
     >                    nlhd,jmtd,nmzd,nmzxyd,natd,nrhomfile,ndomfile,
     >                     neq,invs,invs2,rmt,rmsh,zatom,namat,dx,delz,
     >                     z1,sk3,tau,bmat,ig,ig2,symor,vol,taual,
     >                     jspins,film,nvac,volmts,volint,area,
     >                     ntype,jri,ntypsy,nlh,
     >                     mrot,ng2,ng3,kv3,igfft2,invtab,
     >                     igfft,pgfft2,pgfft,nstr2,nstr,
     >                     kimax,kimax2,nmz,nmzxy,sigma,odi,odl)
               CALL cpu_time(time2)
               CALL outtime(
     >          'gen. spin-up and -down density:',time2-time1)
             ENDIF
c---> pk non-collinear

             CALL cpu_time(time1)
             IF (disp) THEN 
                reap = .false.
             ELSE
                reap = .true.
             ENDIF
             total = .true.
#ifdef CPP_HTML
             CALL html_out('pot',it,itmax,ivers)
#endif
           ENDIF!(pot8)
         ENDIF !irank.eq.0

c--- J<
      IF(l_jenerg) GOTO 234

         alph1(:)=alph(:)
         stop80= .false. 
         IF ( eig66(1) .or. (l_soc .and. l_ss) ) THEN
           IF ( (l_J).or.(nqpt/=1).or.(nmagn/=1).or.(phnd/=1) ) THEN
             STOP 'fleur: J-loop with eig66=T or ss+soc'
           ENDIF
         ENDIF
         DO qcount=1,nqpt
           IF (l_J) THEN
             qss(:)=qj(:,qcount)
             qn = ( qss(1)**2 + qss(2)**2 + qss(3)**2 )
           ENDIF
           IF ( l_J.AND.(irank.EQ.0) ) THEN
             WRITE(6,*) 'qss=(',qss(1),',',qss(2),',',qss(3),')'
             CALL cpu_time(timeq1)
           ENDIF
           IF (.not.pot8) THEN
             CALL cpu_time(time1)
c--- J>
             CALL vgen(
     >             it,reap,vchk,namat,icorr,total,krla,
     >             ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd,
     >             ncvd,nmzd,k1d,k2d,k3d,kxc1d,kxc2d,kxc3d,
     >             nop,n2d,natd,nn2d,nmzxyd,nspd,nq2,nq3,nmzxy,
     >             nmz,mx3,nk2,nk1,ntypsy,lmax,jri,nlh,neq,ncv,nwd,
     >             invs2,invs,zrfs,jspins,ntype,film,nvac,mrot,tau,
     >             z1,area,dvac,delz,omtil,volint,volmts,rmsh,symor,
     >             llh,nmem,mlh,clnu,nsymt,zatom,rmt,dx,pos,amat,
     >             bmat,kv2,kv3,sk2,sk3,ig,ig2,nstr,nstr2,ngopr,
     >             kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,
     >             kimax,igfft,pgfft,kimax2,igfft2,pgfft2,ufft,
     >             pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy,vol,taual,
     >             igrd,ndvgrd,idsprs,isprsv,idsprs0,idsprsl,idsprsi,
     >             idsprsv,chng,sprsv,lwb,nop2,zsigma,sigma,sig_b,
     >             invtab,phi2,odi,ods,odl,odg,plpot,irank,isize,
     <             vr0,te_vcoul,te_veff,te_exc,l_noco,l_ss,qss)

             IF (irank.EQ.0) THEN
             CALL cpu_time(time2)
             CALL outtime(
     >          'generation of potential (total):',time2-time1)

c---> pk non-collinear
c--->          generate the four component matrix potential from spin up
c--->          and down potentials and direction of the magnetic field
               IF (l_noco) THEN
                  CALL cpu_time(time1)
                  CALL vmatgen(
     >                   jspd,k1d,k2d,k3d,n2d,n3d,
     >                   ntypd,ntypsd,nlhd,natd,jmtd,nmzd,nmzxyd,
     >                   ntype,ntypsy,nlh,jri,neq,invs,invs2,
     >                   8,ndomfile,npotmatfile,jspins,film,nvac,
     >                   ng2,ng3,igfft2,
     >                   igfft,pgfft2,pgfft,nstr2,nstr,
     >                   kimax,kimax2,ufft,odi,odl)
                  CALL cpu_time(time2)
                  CALL outtime(
     >             'generation of matrix potential:',time2-time1)
               ENDIF
c---> pk non-collinear
             ENDIF ! irank.eq.0

c
c+t3e
           ENDIF ! .not.pot8
#ifdef CPP_MPI
           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
c
c          ----> eigenvalues and eigenfunctions
c
c--- J<
           IF(l_disp)THEN
             mtypes=1
             nmagn=1
             phnd=1
           ENDIF
         i_J=1
         DO imt=1,mtypes
            DO j_J=i_J,nmagn
             DO phn=1,phnd

               IF (eig66(1)) THEN
                 ! If eig-file exists, use it.
                 ! If eig-file does not exist, create it and stop.
                 !
                 ! eig66(2)=eig66(1)=T <=> in first run of program
                 !                         only 1st variation is done
#ifdef CPP_HDF
                 INQUIRE(file='eig.hdf',exist=eigvar(3))
#else
                 INQUIRE(file=fname(1),exist=eigvar(3))
#endif
                 eigvar(1)= .NOT.eigvar(3)
                 eigvar(2)= eig66(2) .EQV. eigvar(3)
               ELSE
                 eigvar(1)= .true.
                 eigvar(2)= .true.
                 eigvar(3)= .true. 
               ENDIF
               eigvar(2)= eigvar(2) .AND. ( l_soc .AND. (.NOT.l_noco) )
               ! eigvar(1/2)= 1st/2nd var. ; eigvar(3)= calc density,etc  
 
               IF ( l_soc .AND. (.NOT.l_noco) ) THEN
                 evonly(1)= .false.
                 evonly(2)= eonly
               ELSE
                 evonly(1)= eonly
                 evonly(2)= .false.
               ENDIF

               IF ( eigvar(1).OR.eigvar(2) ) THEN
#ifdef CPP_HTML
                 CALL html_out('eig',it,itmax,ivers)
#endif
                 IF (l_J) THEN
                   tkb=0.
                   CALL jcoff(
     >              i_J,j_J,phn,irank,ntypd,natd,ntype,taual,neq,
     >              qss,nmagn,l_magn,l_disp,thetaJ,M,
     X              alph1,alph,beta,l_wr)
                 ENDIF
                 IF (eigvar(1)) THEN 
#ifdef CPP_HPM
                   call f_hpmstart( 5, "eigen" )
#endif
                   ! WRITE(6,fmt='(A)') 'Starting 1st variation ...'
                   CALL eigen(
     >              n3d,memd,nlhd,ntypsd,ntypd,nwdd,natd,
     >              nop,k1d,k2d,k3d,kv2,nkptd,jspd,lmaxd,jmtd,
     >              nvd,nmzxyd,n2d,nmzd,neigd,nv2d,nbasfcn,
     >              irecl,irecl0,irecl1,l_zref,zelec,
     >              it,8,ntype,lpr,form66,evonly(1),secvar,isec1,
     >              film,invs,zrfs,invs2,jspins,nvac,neq,lmax,
     >              bk,wtkpt,rkmax,ellow,elup,nwd,nkpt,
     >              jri,rmt,dx,ntypsy,nlh,nq2,nq3,nmz,nmzxy,
     >              ustep,bbmat,bmat,invsat,invsatnr,invtab,
     >              ig2,ig,rgphs,area,delz,z1,omtil,zatom,namat,
     >              lnonsph,el0,evac0,lepr,irank,isize,
     >              nmem,mlh,llh,rmsh,clnu,ngopr,mrot,taual,
     >              nlod,llod,ello0,llo,nlo,lo1l,l_dulo,ulo_der,lapw_l,
     >              l_noco,l_J,l_ss,alph,beta,qss,npotmatfile,
     >              l_constr,b_con,n_u,lda_u,l_f,fname,l_soc,soc_opt,
     >              gw,sk2,phi2,odi,ods,odd,
     <              e_u_c,nkpt_l,eig_l)
#ifdef CPP_HPM
                   call f_hpmstop( 5 )
#endif
                 ENDIF 
                 IF (eigvar(2)) THEN
                   ! WRITE(6,fmt='(A)') 'Starting 2nd variation ...'
                   CALL eigenso(
     >              jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,nlod,
     >              nwdd,neigd,nvd,lmaxd,nop,natd,ntypsd,
     >              llod,ulo_der,l_dulo,
     >              invsat,invsatnr,lapw_l,bbmat,l_noco,alph,beta,
     >              nbasfcn,film,invs,invs2,nq2,nq3,nvac,jspins,nwd,
     >              irecl1,8,irank,isize,irecl,ntype,soc_opt,nlo,llo,
     >              theta,phi,omtil,neq,lmax,mrot,ngopr,nkpt,nlh,
     >              jri,ntypsy,bmat,rmt,dx,rmsh,taual,invtab,fname,
     >              odi,ods,evonly(2),TYP_REAL)
                 ENDIF 
         
                 CALL cpu_time(time1)
                 IF (irank.EQ.0) THEN
                   IF (it.LT.isec1) THEN
                     CALL outtime('generation of hamiltonian and'//
     >                     ' diagonalization (total):',time1-time2)
                   ELSE
                     CALL outtime('generation of hamiltonian and'//
     >                     ' red. diagonalization (total):',time1-time2)
                   END IF
c+t3e
                 ENDIF
#ifdef CPP_MPI
                 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
                 IF (eig66(1)) THEN
                   IF (irank==0) THEN
                     WRITE (*,fmt='(A)') 
     &                'eig-file created, program stops'
                   ENDIF
                   stop80= .true. 
                 ENDIF

               ENDIF ! ( eigvar(1) .OR. eigvar(2) )

               IF ( eigvar(3) .and. l_soc .and. l_ss ) THEN 
                 CALL ssomat(
     >            jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,nlod,
     >            neigd,nvd,lmaxd,nop,natd,ntypsd,llod,l_dulo,
     >            ulo_der,invsat,invsatnr,lapw_l,bbmat,alph,beta,
     >            nbasfcn,film,invs,invs2,nq2,nq3,nvac,jspins,
     >            nwd,8,irank,isize,irecl,ntype,soc_opt,nlo,llo,
     >            theta,phi,omtil,neq,lmax,mrot,ngopr,nkpt,nlh,
     >            jri,ntypsy,bmat,rmt,dx,rmsh,taual,invtab,
     >            odi,ods,qss,zelec,tkb,sso_opt,TYP_REAL )
                 stop80= .true.
               ENDIF
c-t3e  
c
c              ----> fermi level and occupancies
c

               IF ( eigvar(3) .and. ( .not.(l_soc .and. l_ss) ) ) THEN  
#ifdef CPP_HDF
                 CALL openeig(
     >                 nbasfcn,neigd,nkpt(1),jspd,lmaxd,nlod,ntypd,
     >                 .false.,.true.)
#else
                 IF (l_J) THEN
                   CALL cpu_time(time1)
                   ALLOCATE ( w_iks(neigd,nkptd,jspd) )
                   CALL fermie(
     >                 neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod,
     >                 jspins,nwd,ntype,nkpt,66,irecl,form66,.false., 
     >                 gauss,delgau,zelec,tkb,tria,film,l_noco,l_ss,
     <                 ef,seigscv,ts,w_iks,
     X                 qss,l_J,l_disp,bmat,nkpt_l,eig_l,irank,isize)
                   DEALLOCATE ( w_iks )
                   CALL cpu_time(time2)
                   CALL outtime(
     >              'determination of fermi energy:',time2-time1)
                 ELSE 
                 OPEN (66,file=fname(1),access='direct',
     >                   form='unformatted',recl=irecl,status='old')
                 ENDIF
#endif
                 IF ( l_soc .and. (.not. l_noco) ) neigd = 2*neigd
                 ALLOCATE ( w_iks(neigd,nkptd,jspd) )
c-dw-chs
                 IF ( (irank.EQ.0).AND.(.not.l_J) ) THEN
#ifdef CPP_HTML
                   CALL html_out('fer',it,itmax,ivers)
#endif
                   CALL cpu_time(time1)

                   IF ( l_soc .and. (.not. l_noco) ) THEN
                     DO n=1,nwd
                       zelec(n) = zelec(n)*2
                     ENDDO
                     CALL fermie(
     >                 neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod,
     >                 jspins,nwd,ntype,nkpt,66,irecl,form76,evonly(2),
     >                 gauss,delgau,zelec,tkb,tria,film,l_noco,l_ss,
     <                 ef,seigscv,ts,w_iks,
     X                 qss,l_J,l_disp,bmat,nkpt_l,eig_l,irank,isize)
                     seigscv = seigscv/2
                     ts = ts/2
                     DO n=1,nwd
                       zelec(n) = zelec(n)/2
                     ENDDO
                   ELSE
                     CALL fermie(
     >                 neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod,
     >                 jspins,nwd,ntype,nkpt,66,irecl,form66,.false., 
     >                 gauss,delgau,zelec,tkb,tria,film,l_noco,l_ss,
     <                 ef,seigscv,ts,w_iks,
     X                 qss,l_J,l_disp,bmat,nkpt_l,eig_l,irank,isize)
                   ENDIF
                   CALL cpu_time(time2)
                   CALL outtime(
     >              'determination of fermi energy:',time2-time1)
c+t3e
                 ENDIF

                 IF (eonly) THEN
#ifdef CPP_HDF
                   CALL closeeig()
#else
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
                   IF (l_J) THEN
                   IF (irank==0) CLOSE (66,status='delete')
                   ELSE
                   CLOSE (66,status='delete')
                   ENDIF
#else
                   CLOSE (66)
#endif
#endif
                   IF (.not. l_J) THEN 
                     DEALLOCATE( w_iks )
#ifdef CPP_MPI
                     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
                     IF (irank==0) THEN
                       WRITE (*,fmt='(A)') 
     &                  'eigenvalues written, program stops'
                     ENDIF
                     stop80= .true. 
                   ENDIF
                 ENDIF ! eonly 

               ENDIF ! ( eigvar(3) .and. .not.(l_soc .and. l_ss) ) 
c--- J<
               IF(l_J) THEN
                 IF (.not. eonly) THEN 
                   DEALLOCATE ( w_iks )
                 ENDIF 
                 IF (((i_J.eq.j_J)).OR.(invs.and.(qn.GT.tol))) GOTO 33
               ENDIF
             ENDDO !phn
  33         CONTINUE
            ENDDO !j_J
         i_J=i_J+nmagtype(imt)
         ENDDO !imt
           IF ((irank.EQ.0).AND.(l_J)) THEN
             CALL cpu_time(timeq2)
             CALL outtime('q-point (total)',timeq2-timeq1)
           ENDIF
         ENDDO !qcount
         IF (stop80) THEN
           EXIT ! it
         ENDIF

 234     CONTINUE

         IF (irank.EQ.0) THEN
           IF(l_J) THEN
             IF(.not.l_disp)THEN
               REWIND(113)
               REWIND(114)
               CALL jcoff2(
     >                 natd,taual,nop,mrot,amat,invs,
     >                 ntypd,nmagn,nqpt,nsh,thetaJ,
     >                 mtypes,magtype,nmagtype,neq,film)
             ENDIF
             CLOSE(113)
             CLOSE(114)
           ENDIF
         ENDIF

         IF (.not.l_J) THEN
c--- J>
#ifdef CPP_MPI
         CALL MPI_BCAST(ef,1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(w_iks,neigd*nkptd*jspd,TYP_REAL,0,
     +                                MPI_COMM_WORLD,ierr)
#endif
c-t3e
c
c        ----> initialise force_old
c
         CALL force_0(
     >                 jspd,ntypd,
     X                 force,
     <                 force_old)
c
c        ----> charge density
c
#ifdef CPP_HTML
         CALL html_out('rho',it,itmax,ivers)
#endif 
         IF (dos) THEN
           IF ((it.LT.itmax).AND.(ndir.EQ.-2)) THEN 
             ndir_tmp = 0
             dos_tmp  = .false.
           ELSE
             ndir_tmp = ndir
             dos_tmp  = .true.
           ENDIF 
         ELSE
           ndir_tmp = 0
           dos_tmp = .false.
         ENDIF

!+Wannier stuff for testing purposes
#ifdef CPP_WANN
         l_wann = .false.
         INQUIRE (file='wann_inp',exist=l_wann)
         IF (l_wann) THEN
            CALL wannier(l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >         neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >         nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >         invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >         beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >         jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >         ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >         ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >         z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >         volint,symor,pos,ef,irecl)
            STOP 'wannier done'
         ENDIF
#endif
!-Wannier

         CALL cdngen(
     >    nrhomfile,npotmatfile,noinpfile,kcrel,
     >    cdinf,dos_tmp,frcor,slice,ctail,ndir_tmp,vacdos,integ,layers,
     >    izlay,kk,e1s,e2s,nnne,pallst,layerd,neigd,nkptd,nv2d,nbasfcn,
     >    ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd,
     >    nvd,nmzd,k1d,k2d,k3d,kq1d,kq2d,kq3d,msh,
     >    nop,n2d,natd,nwdd,nmzxyd,nspd,nn3d,nstd,irecl0,
     >    namat,igfft,pgfft,ig,sk3,nstr,nstr2,kv2,enmix,amat,bmat,
     >    mx1,mx2,mx3,film,symor,zrfs,invs2,invs,ngopr,volint,vol,
     >    volmts,zatom,jri,l_soc,lda_u(:)%l,n_u,
     >    jspins,nwd,nkpt,ntype,nq3,nmz,nq2,nvac,nmzxy,
     >    nsymt,clnu,llh,mlh,nlh,nmem,ntypsy,pos,taual,
     >    rmt,rmsh,dx,neq,lmax,z1,omtil,area,delz,theta,phi,
     >    kq1_fft,kq2_fft,kq3_fft,nq3_fft,kmxq2_fft,kmxq_fft,
     >    kv3,nop2,mrot,zelec,tau,ig2,rgphs,irank,isize,
     >    ncst,lchange,lchg_v,w_iks,l_f,l_geo,bbmat,
     >    invarop,invarind,multab,invtab,invsat,invsatnr,d_wgn,
     >    nlod,llod,nlo,llo,l_dulo,ulo_der,
     >    lapw_l,llochg,skiplo,nlol,lo1l,gw,
     >    lepr,sigma,tworkf,nstars,nstm,starcoeff,locx,locy,
     >    fname,sk2,phi2,odi,ods,
     X    l_noco,l_ss,qss,sso_opt,l_mperp,l_relax,l_constr,mix_b,
     X    ello0,el0,evac0,force,alph,beta,b_con,
     <    seigc,e1_dos,e2_dos,ef,sig_dos)
c
         DEALLOCATE ( w_iks )
         IF ( l_soc .and. (.not. l_noco) ) neigd=neigd/2 
c+t3e
#ifdef CPP_MPI
         CALL MPI_BCAST(evac0,2*jspd*nwdd,TYP_REAL,0,
     +                                MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(el0,(lmaxd+1)*ntypd*jspd*nwdd,TYP_REAL,0,
     +                                MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(ello0,nlod*ntypd*jspd,TYP_REAL,0,
     +                                MPI_COMM_WORLD,ierr)
         IF (l_noco) THEN
           DO n= 1,ntype
             IF (l_relax(n)) THEN
               CALL MPI_BCAST(alph(n),1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(beta(n),1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
             ENDIF
           ENDDO
           IF (l_constr) THEN
             CALL MPI_BCAST(b_con,2*ntypd,TYP_REAL,0,
     +                                MPI_COMM_WORLD,ierr)
           ENDIF
         ENDIF
#endif
         IF (irank.EQ.0) THEN
c-t3e
           CALL cpu_time(time1)
           CALL outtime('generation of new charge density (total)',
     +                  time1-time2)
c
           IF (ndir.GT.0) STOP 'NDIR'
c          ----> output potential and potential difference
           IF (disp) THEN
             reap = .false.
             total = .false.
             CALL cpu_time(time2)
             CALL vgen(
     >                it,reap,vchk,namat,icorr,total,krla,
     >                ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd,
     >                ncvd,nmzd,k1d,k2d,k3d,kxc1d,kxc2d,kxc3d,
     >                nop,n2d,natd,nn2d,nmzxyd,nspd,nq2,nq3,nmzxy,
     >                nmz,mx3,nk2,nk1,ntypsy,lmax,jri,nlh,neq,ncv,nwd,
     >                invs2,invs,zrfs,jspins,ntype,film,nvac,mrot,tau,
     >                z1,area,dvac,delz,omtil,volint,volmts,rmsh,symor,
     >                llh,nmem,mlh,clnu,nsymt,zatom,rmt,dx,pos,amat,
     >                bmat,kv2,kv3,sk2,sk3,ig,ig2,nstr,nstr2,ngopr,
     >                kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,
     >                kimax,igfft,pgfft,kimax2,igfft2,pgfft2,ufft,
     >                pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy,vol,taual,
     >                igrd,ndvgrd,idsprs,isprsv,idsprs0,idsprsl,idsprsi,
     >                idsprsv,chng,sprsv,lwb,nop2,zsigma,sigma,sig_b,
     >                invtab,phi2,odi,ods,odl,odg,plpot,0,1, ! irank,isize,
     <                vr0,te_vcoul,te_veff,te_exc,l_noco,l_ss,qss)
             CALL cpu_time(time1)
             CALL outtime(
     +                'generation of potential (total):',time1-time2)
             CALL potdis(
     >                jspd,n3d,n2d,nmzd,nmzxyd,ntypd,nlhd,ntypsd,natd,
     >                jspins,nq3,nq2,nmz,nmzxy,ntype,nlh,ntypsy,
     >                jmtd,jri,dx,rmsh,neq,omtil,area,vol,
     >                film,nvac,zrfs,invs,invs2,delz,k1d,k2d,k3d,
     >                nk1,nk2,nk3,mx3,mx2,mx1,ig2,rgphs,ig,ustep)
           END IF
c
ci         ----> total energy
c
#ifdef CPP_HTML
           CALL html_out('tot',it,itmax,ivers)
#endif
           CALL totale(
     >               ntypd,memd,ntypsd,n3d,jspd,nmzxyd,n2d,
     >               natd,jmtd,nmzd,layerd,nop,nwdd,nlhd,nlod,
     >               it,seigscv,seigc,ts,vr0,te_vcoul,te_veff,te_exc,
     >               invs,invs2,film,nq2,nq3,l_f,nvac,jspins,ntype,
     >               neq,mrot,llh,mlh,nlh,nmem,ngopr,jri,ntypsy,
     >               dx,tau,pos,zatom,rmt,rmsh,clnu,l_noco,e_u_c,n_u,
     >               odi,ods,odd,amat,bmat,
     X               invtab,l_geo,force_old,force)

           CALL cpu_time(time2)
           CALL outtime('determination of total energy:',time2-time1)
           CALL cpu_time(time2)
c+t3e
         ENDIF ! irank.EQ.0 
#ifdef CPP_HDF
         CALL closeeig()
#else
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
         CLOSE (66) ! ,status='delete')
#else
         CLOSE (66)
#endif
#endif
         ENDIF !(if not l_J)
         ELSE
            alpha = alpha - 10.
         END IF !(if alpha <10.)
         IF (.NOT.l_J) THEN

         IF (irank.EQ.0) THEN
c-t3e
c
c          ----> mix input and output densities
c
#ifdef CPP_HTML
           CALL html_out('mix',it,itmax,ivers)
#endif

           CALL mix(
     >            k1d,k2d,k3d,ntypd,ntypsd,jmtd,jspd,nmzd,nmzxyd,n3d,
     >            nrhomfile,film,invs2,invs,nvac,jspins,nq2,nq3,nmzxy,
     >            nmz,z1,delz,ntype,rmt,zatom,dx,rmsh,neq,ntypsy,jri,
     >            nlh,imix,maxiter,alpha,spinf,it,namat,nlhd,n2d,
     >            natd,nop,kv3,nstr,ig,ig2,bmat,sk3,mrot,tau,invtab,
     >            symor,volint,taual,volmts,l_noco,alph,beta,
     >            kimax,igfft,pgfft,ufft,nstr2,area,omtil,sigma,
     >            vol,odi)
c
           CALL cpu_time(time1)
           CALL outtime('mixing:',time1-time2)
           WRITE (6,FMT=8130) it
           WRITE (16,FMT=8130) it
 8130      FORMAT (/,5x,'******* it=',i3,'  is completed********',/,/)
           CALL outtime('ITERATION (total):',time1-start_it)
c+t3e
         ENDIF ! irank.EQ.0
c--- J<
         ELSE
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
         CLOSE (66) ! ,status='delete')
#else
         CLOSE (66)
#endif
         ENDIF !(if not l_J)
c--- J>

#ifdef CPP_MPI
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
c-t3e
c+fo
         l_endit=.false.
         INQUIRE (file='inp_new',exist=l_endit)
         IF (l_endit) STOP ' GEO new inp created ! '
c-fo
#if (defined(DUMP_DATA_OLD) || defined(DUMP_DATA_HDF))
         print*, 'SINGLE_ITER: stopping after first iteration'
         exit
#endif
   80 CONTINUE
      IF (irank == 0 ) THEN
        DEALLOCATE ( pgft2x,pgft2y,pgft2xx,pgft2xy,pgft2yy,
     +               pgft1x,pgft1y,pgft1xx,pgft1xy,pgft1yy) 
      ENDIF
#ifdef DUMP_DATA_HDF
! close hdf5, probably conflicts with CPP_HDF, which is not used
      call binmat_end()
#endif
c+t3e
#ifdef CPP_HDF
#ifdef CPP_HPM
      call f_hpmstop( 1 ) 
      call f_hpmterminate( taskID ) 
#endif
      CALL hdf_close()
#endif
#ifdef CPP_MPI
      CALL MPI_TYPE_FREE(TYP_REAL,ierr)
      CALL MPI_TYPE_FREE(TYP_COMPLEX,ierr)
      CALL MPI_FINALIZE(ierr)
#endif
#ifdef CPP_HTML
      CALL html_out('bot',0,0,ivers)
#endif
c-t3e
!      IF (l_J) 
      DEALLOCATE (eig_l)
      DEALLOCATE ( clnu,ustep,ig,ig2,igz,jri,kv2,kv3,llh,lmax,mlh,mrot,
     +             ncv,neq,ngopr,nlh,nmem,nstr2,ntypsy,nstr,
     +             igfft,igfft2,nflip,kq1_fft,kq2_fft,kq3_fft,nq2_fft,
     +             nq3_fft,kmxq2_fft,kmxq_fft,ncst,izlay,invarop,
     +             invarind,multab,invtab,invsat,invsatnr,lnonsph,nkpt,
     +             dx,pos,tau,rmsh,rmt,sk2,sk3,taual,volmts,zatom,el0,
     +             evac0,rgphs,force,force_old,bk,wtkpt,ellow,elup,
     +             rkmax,zelec,pgfft,pgfft2,ufft,bmu,vr0,
     +             lchange,lchg_v,l_geo,nlo,
     +             llo,skiplo,ello0,llochg,lo1l,nlol,alph,beta,M,l_magn,
     +             l_relax,lda_u,lapw_l,b_con,l_dulo,enmix,d_wgn,
     +             ulo_der,soc_opt,ig1,kv1,nstr1,ngopr1,mrot1,tau1,
     +             magtype,nmagtype,invtab1,multab1,igfft1,pgfft1)

      STOP 'all done'
      END PROGRAM fleur
