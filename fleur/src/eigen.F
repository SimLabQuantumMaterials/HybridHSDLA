      MODULE m_eigen
      CONTAINS
      SUBROUTINE eigen(
     >                 n3d,memd,nlhd,ntypsd,ntypd,nwdd,natd,
     >                 nop,k1d,k2d,k3d,kv2,nkptd,jspd,lmaxd,jmtd,
     >                 nvd,nmzxyd,n2d,nmzd,neigd,nv2d,nbasfcn,
     >                 irecl,irecl0,irecl1,l_zref,zelec,
     >                 it,nu,ntype,lpr,form66,eonly,secvar,isec1,
     >                 film,invs,zrfs,invs2,jspins,nvac,neq,lmax,
     >                 bk,wtkpt,rkmax,ellow,elup,nwd,nkpt,
     >                 jri,rmt,dx,ntypsy,nlh,nq2,nq3,nmz,nmzxy,
     >                 ustep,bbmat,bmat,invsat,invsatnr,invtab,
     >                 ig2,ig,rgphs,area,delz,z1,omtil,zatom,namat,
     >                 lnonsph,el0,evac0,lepr,irank,isize,
     >                 nmem,mlh,llh,rmsh,clnu,ngopr,mrot,taual,nlod,
     >                 llod,ello0,llo,nlo,lo1l,l_dulo,ulo_der,lapw_l,
     >                 l_noco,l_J,l_ss,alph,beta,qss,npotmatfile,
     >                 l_constr,b_con,n_u,lda_u,l_f,fname,l_soc,soc_opt,
     >                 gw,sk2,phi2,odi,ods,odd,
     <                 e_ldau,nkpt_l,eig_l)
c
c*********************************************************************
c     sets up and solves the eigenvalue problem for a basis of lapws.
c     as written, each k-point could be done as a separate task in a
c     multi-tasking environment with simple and obvious changes.
c     it = iteration number in this run
c     nu = file number of potential
c     ntype = number of distinct atom types
c     lpr > 0 prints eigenvectors also
c     form66 = .true. for formatted eigenvector file
c     eonly =.true. calculates eigenvalues only 
c                  m. weinert  jan. 1987
c
c nv,   nvd     ... actual length & dimension of EV without LO's 
c nmat, nbasfcn                                   including LO's
c***********************************************************************
c     Definition of irecl1:
c
c  (1)l_zref=.true.
c     REAL:   bkpt(3),rk(nvd)                              ->(nvd+3)  *8
c     INTEGER:k1(nvd),k2(nvd),k3(nvd),nv,nred,matind(nvd,2)->(5*nvd+2)*4
c                                                 (oder *8 auf der cray)
c  (2)l_zref=.false.
c     without nred,matind(nvd,2) -> -(2*nvd+1)*4 (or 8 on the  cray)
c
c               Stefan Bl"ugel and Gustav Bihlmayer, IFF, FZJ, 31.Dec.97
c**********************************************************************
c     modified for explicit use of z-reflection symmetry in franza 
c---------------------------------------------------
c this routine is called by: fleur.F
c
c eigen. --+-- loddop       ... read potential
c          +-- lodpot
c          |
c         +-+               ... window loop
c         | +-- tlmplm
c         |+-+              ... k-point loop
c         || +-- apws
c         || +-- hsint
c         || +-- (hssphn)
c         || +-- hsvac - vacfun - vacu(d)z
c         || +-- franza || aline
c         || +-- (outeig_muff)
c         || +-- (aline_muff)
c         || +-- outeig
c         |+-+
c         +-+
c---------------------------------------------------
c        g. bihlmayer '96                                       
c**********************************************************************
      USE m_constants, ONLY : pimach
      USE m_lodpot
      USE m_tlmplm
      USE m_apws
      USE m_hssphn
      USE m_zsymsecloc
      USE m_aline
      USE m_alinemuff
      USE m_hsint
      USE m_hsvac
      USE m_od_hsvac
      USE m_usetup
      USE m_types,    ONLY : t_utype
      USE m_loddop
      USE m_outeig
      USE m_od_types, ONLY : od_inp, od_sym, od_dim
      USE m_cputime
      USE m_outtime
#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : openeig, writeeig, closeeig, writebas
#endif
#ifdef CPP_MPI
      USE m_mpimakegroups
#endif
      USE m_gweig
      USE binmat ! for matrix dumping

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT(IN) :: it,lpr,ntype,nu,nwd,isec1,irank,isize
      INTEGER,INTENT(IN) :: n3d,memd,nlhd,ntypsd,ntypd,nwdd,natd
      INTEGER,INTENT(IN) :: nop,k1d,k2d,k3d,nkptd,jspd,lmaxd,jmtd
      INTEGER,INTENT(IN) :: nvd,nmzxyd,n2d,nmzd,neigd,nv2d,nkpt_l
      INTEGER,INTENT(IN) :: nq2,nq3,nmz,nmzxy,lepr,jspins,nvac,irecl
      INTEGER,INTENT(IN) :: nlod,llod,npotmatfile,nbasfcn,n_u,gw
      REAL,   INTENT(IN) :: omtil,area,delz,z1
      LOGICAL,INTENT(IN) :: eonly,secvar,film,invs,zrfs,invs2
      LOGICAL,INTENT(IN) :: l_noco,l_J,l_ss,l_constr,l_f,l_soc
      LOGICAL,INTENT(INOUT) :: form66
      REAL,   INTENT(OUT):: e_ldau
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: ustep(n3d)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: lnonsph(ntypd),nkpt(nwdd),lmax(ntypd)
      INTEGER, INTENT (IN) :: invsat(natd),ngopr(natd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop),invsatnr(natd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd),kv2(2,n2d)
      INTEGER, INTENT (IN) :: lo1l(0:llod,ntypd),lapw_l(ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: ello0(nlod,ntypd,jspd),zatom(ntypd)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: bk(3,nkptd,nwdd),wtkpt(nkptd,nwdd)
      REAL,    INTENT (IN) :: ellow(nwdd),elup(nwdd),rkmax(nwdd)
      REAL,    INTENT (IN) :: rmt(ntypd),dx(ntypd),evac0(2,jspd,nwdd)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3),taual(3,natd)
      REAL,    INTENT (IN) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd),qss(3)
      REAL,    INTENT (IN) :: b_con(2,ntypd),zelec(nwdd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      LOGICAL, INTENT (IN) :: soc_opt(ntypd+2)
      CHARACTER*2, INTENT    (IN) :: namat(0:103)
      CHARACTER*12,INTENT    (IN) :: fname(3)
      TYPE (t_utype), INTENT (IN) :: lda_u(ntypd)
c-odim
      REAL,    INTENT (IN) :: sk2(n2d),phi2(n2d)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
      REAL,   INTENT(OUT):: eig_l(neigd+5,nkpt_l)
C     ..
C     .. Local Scalars ..
      INTEGER irecl0,irecl1,jsp,nb0,nk,nk0,nmat,nw,nred,ndim
      INTEGER iter,ne,nblw,matsize,lmd,lmplmd,nrec,matel,lh0
      INTEGER loplod,nspins,isp,n_u_l,nlotot,l,i,j,ne_sv,err,gwc
      INTEGER mlotot,mlolotot,mlot_d,mlolot_d,nlot_d
      LOGICAL l_zref,l_wu
      REAL cp_aline,cp_apws,cp_axbxt,cp_hsint,cp_hssphn,cp_hsvac,
     +     cp_outeig,cp_outeigm,cp_tlmplm,cp_seclr4,time1,time2,
     +     cp_sph,cp_nsph,sfp
      REAL evac_sv(jspd)
C     .. for matrix dumping ..
      character(len=32) :: hamilton_name = 'xxx_Hmat.bin', 
     +                      overlap_name = 'xxx_Smat.bin',
     +                      misc_name
C     ..
C     .. Local Arrays ..
      INTEGER, PARAMETER :: lmaxb=3
      INTEGER, ALLOCATABLE :: matind(:,:),kveclo(:)
      INTEGER, ALLOCATABLE :: nv(:),nv2(:),k1(:,:),k2(:,:),k3(:,:)
      REAL,    ALLOCATABLE :: bkpt(:),ddn(:,:,:),uloulopn(:,:,:,:)
      REAL,    ALLOCATABLE :: eig(:),rk(:,:),ello(:,:,:),evac(:,:,:)
      REAL,    ALLOCATABLE :: el(:,:,:,:),bound_lo(:),bound_up(:)
      REAL,    ALLOCATABLE :: us(:,:,:),uds(:,:,:),dus(:,:,:)
      REAL,    ALLOCATABLE :: duds(:,:,:),ulos(:,:,:),dulos(:,:,:)
      REAL,    ALLOCATABLE :: uulon(:,:,:),dulon(:,:,:)
#ifdef CPP_INVERSION
      REAL,    ALLOCATABLE :: z(:,:),a(:),b(:)
#else
      COMPLEX, ALLOCATABLE :: z(:,:),a(:),b(:)
#endif
      COMPLEX, ALLOCATABLE :: vs_mmp(:,:,:,:),vpw(:,:),vzxy(:,:,:,:)
      COMPLEX, ALLOCATABLE :: vpwtot(:,:)
      REAL,    ALLOCATABLE :: vz(:,:,:),vr(:,:,:,:),vr0(:,:,:)
      REAL,    ALLOCATABLE :: vrtot(:,:,:,:)
      COMPLEX, ALLOCATABLE :: tuu(:,:,:),tdd(:,:,:),tud(:,:,:)
      COMPLEX, ALLOCATABLE :: tuulo(:,:,:,:),tdulo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: tdu(:,:,:),tuloulo(:,:,:,:)

      INTEGER, ALLOCATABLE :: ind(:,:,:,:)
      
      CHARACTER*8 dop,iop,name(10)
C     ..
C     .. External Subroutines ..
      EXTERNAL outeig_muff
c
      INTEGER n_start,n_groups,n_rank,n_size,n,n_stride
      INTEGER SUB_COMM
#if (defined(CPP_MPI) && !defined(CPP_T90))
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif
#ifdef DUMP_DATA_OLD
      integer :: dump_i
#endif
#ifdef DUMP_DATA_MT_HSSPHN_ONLY
      logical :: dump_hs
#endif
!
! --> Allocate for the small stacks:
!
      ALLOCATE ( uloulopn(nlod,nlod,ntypd,jspd),nv(jspd),nv2(jspd) )
      ALLOCATE ( ddn(0:lmaxd,ntypd,jspd),eig(neigd),bkpt(3) )
      ALLOCATE ( us(0:lmaxd,ntypd,jspd),uds(0:lmaxd,ntypd,jspd) )
      ALLOCATE ( dus(0:lmaxd,ntypd,jspd),duds(0:lmaxd,ntypd,jspd) )
      ALLOCATE ( ulos(nlod,ntypd,jspd),dulos(nlod,ntypd,jspd) )
      ALLOCATE ( uulon(nlod,ntypd,jspd),dulon(nlod,ntypd,jspd) )
      ALLOCATE ( bound_lo(nwdd),bound_up(nwdd),ello(nlod,ntypd,jspd) )
      ALLOCATE ( el(0:lmaxd,ntypd,jspd,nwdd),evac(2,jspd,nwdd) )
      ALLOCATE ( k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd),rk(nvd,jspd) )
c
c --> some parameters first
c
      lmd     = lmaxd* (lmaxd+2)
      lmplmd  = (lmd* (lmd+3))/2
      loplod  = nlod * (nlod+1)/2
c
c     determine the total number of lo's : nlotot 
c
      nlotot = 0 ; mlotot = 0 ; mlolotot = 0
      DO n = 1, ntype
        mlotot = mlotot + nlo(n)
        mlolotot = mlolotot + nlo(n)*(nlo(n)+1)/2
        DO l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        ENDDO
      ENDDO
      nlot_d = max(nlotot,1)
      ALLOCATE ( kveclo(nlot_d) )
C     ..
c---> timing
      IF (irank.eq.0) WRITE (2,8005)
 8005 FORMAT ('EIGENVALUE PART (eigen):')
      cp_aline = 0.0
      cp_apws = 0.0
      cp_axbxt = 0.0
      cp_hsint = 0.0
      cp_hssphn = 0.0
      cp_sph = 0.0
      cp_nsph = 0.0
      cp_hsvac = 0.0
      cp_outeig = 0.0
      cp_outeigm = 0.0
      cp_tlmplm = 0.0
      cp_seclr4 = 0.0
!
!  if gw = 1 or 2, we are in the first or second run of a GW  calculation
!  if gw = 1 we just proceed as normal (one round),
!  if gw = 2 it's the second run: write out the eigenfunctions and
!  the matrix elements with the xc-potential (needs two rounds)
!
      gwc = 1
      IF (form66) THEN
         IF (irank.eq.0) 
     +   OPEN (66,file='eig',form='formatted',status='unknown')
      ELSE
#ifdef CPP_HDF
         CALL openeig(
     >                nbasfcn,neigd,nkpt(1),jspd,lmaxd,nlod,ntypd,
     >                .true.,.false.)
#else
         IF (l_J) THEN
         IF (irank.eq.0)
     +   OPEN (66,file=fname(1),access='direct',form='unformatted',
     +                              recl=irecl,status='unknown')
         ELSE
         OPEN (66,file=fname(1),access='direct',form='unformatted',
     +                              recl=irecl,status='unknown')
         ENDIF
#endif
      ENDIF
c
c look, if WU diagonalisation
c
      IF (it.LT.isec1) THEN
         IF (irank.eq.0) WRITE (6,FMT=8110) it,isec1
 8110    FORMAT (' IT=',i4,'  ISEC1=',i4,' standard diagonalization')
         l_wu = .false.
      ELSE
         IF (irank.eq.0) WRITE (6,FMT=8120) it,isec1
 8120    FORMAT (' IT=',i4,'  ISEC1=',i4,' reduced diagonalization')
         l_wu = .true.
         IF (form66) STOP 'Wu.and.form66'
      END IF
c
c load potential from file pottot (=unit 8)
c
      ALLOCATE ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
      ALLOCATE ( vz(nmzd,2,4), vr(jmtd,0:nlhd,ntypd,jspd) )
      ALLOCATE ( vr0(jmtd,ntypd,jspd) )
      OPEN (nu,file='pottot',form='unformatted',status='old')
      IF (gw.eq.2) THEN
        ALLOCATE ( vpwtot(n3d,jspd), vrtot(jmtd,0:nlhd,ntypd,jspd) )
        WRITE(6,'(A/A/A/A)')
     &  'Info: vxc matrix elements for GW will be calculated in gw_vxc',
     &  'Info: and stored in "vxc", the values obtained from the',
     &  'Info: original implementation are saved to "vxc.old".'
      ENDIF
 999  CONTINUE
      CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nu,natd,neq,
     <            iop,dop,iter,vr,vpw,vz,vzxy,name)
      CLOSE(nu)
      
!
! some modifications for gw-calculations
!
      IF (gw.eq.2.and.gwc.eq.1) THEN
        vrtot  = vr  ! store potential for subroutine gw_vxc
        vpwtot = vpw ! 
      ENDIF

      IF (gwc==1) THEN
         vr0(:,:,:) = vr(:,0,:,:) 
         lh0 = 1
      ELSE IF (gwc==2) THEN
         lh0 = 0                         ! for a gw-calculation, we
         sfp = sqrt(4.0*pimach())        ! now evaluate matrix elements
         DO jsp = 1,jspins               ! with the coulomb potential
           DO n = 1,ntype                ! but with explicit kinetic energy
             DO j = 1,jri(n)
               vr(j,0,n,jsp) = vr(j,0,n,jsp)-vr0(j,n,jsp)*sfp/rmsh(j,n)
             ENDDO
           ENDDO
         ENDDO
      ENDIF
c
c set energy parameters (normally to that, what we read in)
c
      IF (gwc /= 2) THEN
      CALL lodpot(
     >            irank,jspd,lmaxd,jmtd,nlhd,ntypd,nwdd,nmzd,
     >            jspins,ntype,nwd,film,nvac,nmz,lepr,
     >            jri,dx,rmt,rmsh,lmax,vr,vz,llo,zatom,
     >            el0,evac0,ello0,nlo,nlod,l_dulo,ellow,elup,
     <            el,evac,ello,bound_lo,bound_up)
      ENDIF
#ifdef DUMP_DATA_OLD
! save energy parameters
      !if (ntypd > 1) stop "Dumping only works for a single atom type"
      open(321, file='scalars.dat', status='old', position='append')
      write(321, *) el(:, 1, 1, 1)
      close(321)
#endif
#ifdef DUMP_DATA_HDF
! HDF dump
! * energy parameters for spin direction 1 and window 1
      call write_hdf_mat('energy_params', el(:, :, 1, 1))
! * potential (spherically symmetric (?) part for radial solution)
!   for spin direction 1
      call write_hdf_mat('potential', vr0(:, :, 1),
     +"Rectangular shape, actual elements need to be "//
     +"extracted per-type using the number of "//
     +"mesh points; dimensions are (types, max_mesh_pts)")
#endif
c
c---> open direct access files for the t(l'm',lm) matrices (unit 28)
c---> and the lapw lists (unit 29) (note a cray word has 8 bytes)
c
      IF (l_f.AND.(irank.EQ.0)) THEN
      OPEN (28,file=fname(2),access='direct',recl=irecl0)
      IF (jspins.EQ.2) THEN
        OPEN (38,file=fname(3),status='unknown',access='direct',
     +        recl=irecl0)
      ENDIF
      ENDIF

      nb0 = 0
      nk0 = 0

c---> set up and solve the eigenvalue problem
c---> loop over energy windows
#ifdef DISPLAY_INFO
! display information about iteration parameters
      print "('#', I2, ': ', '#windows = ', I, ', #spins = ', I)", it, nwd, nspins
      print*, 'nkpt: ', nkpt, ', ntypd: ', ntypd
#endif

      DO 30 nw = 1,nwd
c
#if ( defined(CPP_MPI) && defined(CPP_EVP) )
c
         CALL mpi_make_groups(
     >                        irank,isize,nbasfcn,nkpt(nw),neigd,
     >                        zelec(nw),nw,nwd,ntypd,nvd,lmaxd,lmd,
     >                        nlod,llod,mlotot,mlolotot,jspd,n_u,
     >                        lmplmd,l_noco,l_ss,l_soc,l_constr,
     <                        n_start,n_groups,n,matsize,ne,
     <                        n_rank,n_size,SUB_COMM)
         n_stride = nkpt(nw)/n_groups
         IF (n_size > 1) l_zref=.false.
c
#else
c
         n_rank = 0
         n_size = 1
#ifdef CPP_MPI
         n_start = irank + 1
         n_stride = isize
#else
         n_start = 1
         n_stride = 1
#endif
         n = nbasfcn/n_size
         matsize = nbasfcn * (nbasfcn+1)/2
         ne = neigd
#endif
         ne_sv = ne
         ALLOCATE ( a(matsize), stat = err )
         IF (err.NE.0) THEN
           WRITE (*,*) 'eigen: an error occured during allocation of'
           WRITE (*,*) 'the Hamilton Matrix: ',err,'  size: ',matsize
           STOP
         ENDIF
         ALLOCATE ( b(matsize), stat = err )
         IF (err.NE.0) THEN
           WRITE (*,*) 'eigen: an error occured during allocation of'
           WRITE (*,*) 'the overlap Matrix: ',err,'  size: ',matsize
           STOP
         ENDIF
         ALLOCATE (  matind(nbasfcn,2) )
c
c--->    loop over spins
         nspins = jspins
         IF (l_noco) nspins = 1
C
C        Append information about file eig to gwa
         IF(gw.eq.2.and.gwc.eq.1) THEN
           IF(nwd.gt.1) STOP 'eigen: Only use one energy window for GW!'
           OPEN(15,file='gwa',status='old',form='unformatted')
           READ(15)
           READ(15)
           READ(15)
           WRITE(15)
     &       irecl,n_start,n_stride,n_rank,n_size,nvd,nbasfcn,nlotot
           CLOSE(15)
         ENDIF
C  ..
C  LDA+U
         IF ((n_u.GT.0).AND.(nw.EQ.nwd)) THEN
            ALLOCATE( vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins) )
            n_u_l = n_u
            CALL u_setup(
     >              jspd,ntype,lmaxd,lmaxb,jmtd,nlhd,nw,n_u_l,jspins,
     >                   lda_u,neq,el(0,1,1,nw),vr,jri,dx,rmsh,irank,
     <                   vs_mmp,e_ldau)
         ELSE
            ALLOCATE( vs_mmp(-lmaxb:-lmaxb,-lmaxb:-lmaxb,1,2) )
         ENDIF
c
c--->    loop over k-points: each can be a separate task
         DO 20 jsp = 1,nspins
c+do
            nb0=(nw-1)*ntype
c-do
c--->       set up k-point independent t(l'm',lm) matrices
c
            CALL cpu_time(time1) ; err=0
            j = 1 ; IF (l_noco) j = 2
            ALLOCATE(tuu(0:lmplmd,ntypd,j),stat=err)
            ALLOCATE(tud(0:lmplmd,ntypd,j),stat=err)
            ALLOCATE(tdd(0:lmplmd,ntypd,j),stat=err)
            ALLOCATE(tdu(0:lmplmd,ntypd,j),stat=err)
            mlot_d = max(mlotot,1) ; mlolot_d = max(mlolotot,1)
            ALLOCATE(tdulo(0:lmd,-llod:llod,mlot_d,j),stat=err)
            ALLOCATE(tuulo(0:lmd,-llod:llod,mlot_d,j),stat=err)
            ALLOCATE(tuloulo(-llod:llod,-llod:llod,mlolot_d,j),stat=err)
            ALLOCATE(ind(0:lmd,0:lmd,ntypd,j),stat=err )
            IF (err.NE.0) THEN
              WRITE (*,*) 'eigen: an error occured during allocation of'
              WRITE (*,*) 'the tuu, tdd etc.s: ',err,'  size: ',mlotot 
              STOP
            ENDIF
            CALL tlmplm(
     >                  memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd,
     >                  ntype,dx,rmsh,jri,lmax,ntypsy,nwdd,natd,
     >                  lnonsph,lmd,lmplmd,clnu,mlh,nmem,llh,nlh,neq,
     >                  jsp,nw,nb0,irank,mlotot,mlolotot,
     >                  vr(1,0,1,jsp),vr0(1,1,jsp),gwc,lh0,el,secvar,
     >                  nlod,llod,loplod,ello,llo,nlo,lo1l,l_dulo,
     >                  ulo_der,l_f,lmaxb,n_u,lda_u(:)%l,
     >                  vs_mmp(-lmaxb,-lmaxb,1,jsp),
     <                  tuu(0,1,1),tud(0,1,1),tdu(0,1,1),
     <                  tdd(0,1,1),ind(0,0,1,1),tuulo(0,-llod,1,1),
     <                  tdulo(0,-llod,1,1),tuloulo(-llod,-llod,1,1),
     <                  us(0,1,jsp),dus(0,1,jsp),
     <                  uds(0,1,jsp),duds(0,1,jsp),ddn(0,1,jsp),
     <                  ulos(1,1,jsp),dulos(1,1,jsp),uulon(1,1,jsp),
     <                  dulon(1,1,jsp),uloulopn(1,1,1,jsp))
            CALL cpu_time(time2)
            cp_tlmplm = cp_tlmplm + time2 - time1
#ifdef DUMP_DATA_HDF
            call write_hdf_tmat(tuu, tdd, tud, tdu, ind, lmax)
! * hdf: save k-meshes: calculate complete g-mesh for all k-points (cheap)
            call prepare_hdf_kpts(nkpt(1), nvd)
            do nk = 1, nkpt(1) ! window 1
                CALL apws(
     >                   jspd,nvd,nbasfcn,jspins,l_ss,qss,
     >                   bk(1,nk,nw),rkmax(nw),nk,nk0,bbmat,bmat,l_zref,
     >                   n_size,jsp,
     <                   bkpt,nv,rk,k1,k2,k3,matind,nred)
                !foo
                call add_hdf_kpt(nk, bkpt, nv(1),
     +                           k1(:,1), k2(:,1), k3(:,1))
            end do
            call finish_hdf_kpts()
#endif
#ifdef DUMP_DATA_OLD
! * dump t-matrices
            write(misc_name, '(A, ".bin")') 'tmat_uu'
            call write_tlm(tuu, ind, lmax, misc_name)
            write(misc_name, '(A, ".bin")') 'tmat_ud'
            call write_tlm(tud, ind, lmax, misc_name)
            write(misc_name, '(A, ".bin")') 'tmat_du'
            call write_tlm(tdu, ind, lmax, misc_name)
            write(misc_name, '(A, ".bin")') 'tmat_dd'
            call write_tlm(tdd, ind, lmax, misc_name)
#endif

c---> pk non-collinear
c--->       call tlmplm again for the second spin direction in
c--->       each MT, because the t-matrices are needed for both
c--->       spins at once in hssphn
            IF (l_noco) THEN
               isp = 2
               CALL cpu_time(time1)
               CALL tlmplm(
     >                  memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd,
     >                  ntype,dx,rmsh,jri,lmax,ntypsy,nwdd,natd,
     >                  lnonsph,lmd,lmplmd,clnu,mlh,nmem,llh,nlh,neq,
     >                  isp,nw,nb0,irank,mlotot,mlolotot,
     >                  vr(1,0,1,isp),vr0(1,1,isp),gwc,lh0,el,secvar,
     >                  nlod,llod,loplod,ello,llo,nlo,lo1l,l_dulo,
     >                  ulo_der,l_f,lmaxb,n_u,lda_u(:)%l,
     >                  vs_mmp(-lmaxb,-lmaxb,1,isp),
     <                  tuu(0,1,isp),tud(0,1,isp),tdu(0,1,isp),
     <                 tdd(0,1,isp),ind(0,0,1,isp),tuulo(0,-llod,1,isp),
     <                  tdulo(0,-llod,1,isp),tuloulo(-llod,-llod,1,isp),
     <                  us(0,1,isp),dus(0,1,isp),
     <                  uds(0,1,isp),duds(0,1,isp),ddn(0,1,isp),
     <                  ulos(1,1,isp),dulos(1,1,isp),uulon(1,1,isp),
     <                  dulon(1,1,isp),uloulopn(1,1,1,isp))
               CALL cpu_time(time2)
               cp_tlmplm = cp_tlmplm + time2 - time1
            ENDIF
c
            DO 10 nk = n_start,nkpt(nw),n_stride
               nrec = nspins*nk0 + nkpt(nw)*(jsp-1) + nk
               nrec = n_size*(nrec-1) + n_rank + 1
c
c--->         set up lapw list
c
               CALL cpu_time(time1)
               CALL apws(
     >                   jspd,nvd,nbasfcn,jspins,l_ss,qss,
     >                   bk(1,nk,nw),rkmax(nw),nk,nk0,bbmat,bmat,l_zref,
     >                   n_size,jsp,
     <                   bkpt,nv,rk,k1,k2,k3,matind,nred)
               CALL cpu_time(time2)
               cp_apws = cp_apws + time2 - time1
c
               IF (l_noco) THEN
c--->          the file potmat contains the 2x2 matrix-potential in
c--->          the interstitial region and the vacuum
                  OPEN (npotmatfile,FILE='potmat',FORM='unformatted',
     +                  STATUS='old')
               ENDIF
c
c--->         set up interstitial hamiltonian and overlap matrices
c
               CALL cpu_time(time1)
               CALL hsint(
     >              npotmatfile,l_noco,l_J,matsize,jspd,nvd,n3d,k1d,k2d,
     >              k3d,nq3,ig,rgphs,ustep,vpw(1,jsp),k1,k2,k3,nv,jsp,
     >              n_size,n_rank,rk,bk(1,nk,nw),bbmat,nlotot,
     <              nmat,a,b)
               CALL cpu_time(time2)
               cp_hsint = cp_hsint + time2 - time1
#ifdef DUMP_DATA_OLD
! Hamiltonian (a) and Overlap (b) after the interstitial part
               write(misc_name, '(A, ".bin")') 'Hmat_INT'
               call write_mat(a, misc_name)
               write(misc_name, '(A, ".bin")') 'Smat_INT'
               call write_mat(b, misc_name)
#endif               
c
c--->         update with sphere terms
c
#ifdef DUMP_DATA_OLD
! dump parameters for HSSPHN to verify against custom implementation
! we're inside the k-point loop. each kpt has different 
! k1,2,3, (actually G), length stored in nv. 
! additionally, there might be spins, but that's not handled (yet)
              if (nspins > 1) stop 'Data dumping: only 1 spin allowed'
              print*, 'nk=', nk, ', nv(1) = ', nv(1), ', kpt: ', bkpt
! the k-points themselves can be read from kpts (file), as well as the scale
! factor in the 1st line. thus, store the g-vectors from k1, k2, k3 with the
! index of the corresponding k-point, which is 'nk'
              write(misc_name, '(A, "_", I3.3, ".bin")') 'gpts', nk
              call write_mat(k1(:nv(1),1), k2(:nv(1),1), k3(:nv(1),1),
     +                       misc_name)
! write out the rotation matrices as a temporary fix to get the phases right
              do dump_i = 1, natd
                write(misc_name, '(A, "_", I3.3, ".bin")') 'rotmat',
     +                   dump_i
                call write_mat(real(mrot(:, :, invtab(ngopr(dump_i)))),
     +                              misc_name)
              end do
#endif
#ifdef DUMP_DATA_MT_HSSPHN_ONLY
! debug: set H, S to 0, so they only contain the hssphn terms
      a = 0
      b = 0
#endif
               IF (.not.l_wu) THEN
                  CALL cpu_time(time1)
                  CALL hssphn(
     >              nvd,lmaxd,ntypd,nlhd,jmtd,natd,nop,jspd,nwdd,
     >              llod,nlod,loplod,matsize,lmd,lmplmd,nlotot,
     >              ntype,neq,lmax,jri,rmt,rmsh,dx,vr,taual,
     >              bmat,mrot,ngopr,omtil,
     >              nb0,jsp,jspins,nw,secvar,el,ello,lnonsph,
     >              us,dus,uds,duds,ddn,nlo,llo,l_dulo,ulo_der,
     >              kveclo,ulos,dulos,uulon,dulon,uloulopn,lapw_l,
     >              bkpt,k1,k2,k3,rk,nv,invsat,invtab,
     >              l_noco,alph,beta,l_constr,b_con,l_ss,qss,
     >              l_soc,soc_opt,n_u_l,vs_mmp,lmaxb,lda_u(:)%l,odi,ods,
     >              tuu,tud,tdu,tdd,ind,tuulo,tdulo,tuloulo,
     >              mlotot,mlolotot,n_size,n_rank,SUB_COMM,gwc,irank,
     X              a,b,nmat,cp_sph,cp_nsph)

                  CALL cpu_time(time2)
                  cp_hssphn = cp_hssphn + time2 - time1
                  print "('Finished hssphn, cp_hssphn = ', I2 , 
     +             ':', I2, ':', F5.2, '; (', F15.3, 's)')",
     +          int(cp_hssphn)/3600, mod(int(cp_hssphn), 3600)/60,
     +          mod(cp_hssphn, 60.), cp_hssphn
#ifdef DUMP_DATA_OLD
! H and S after _adding_ the MT contribution
                  write(misc_name, '(A, ".bin")') 'Hmat_MT'
                  call write_mat(a, misc_name)
                  write(misc_name, '(A, ".bin")') 'Smat_MT'
                  call write_mat(b, misc_name)
#endif      
#if (defined(DUMP_DATA_HDF) && defined(DUMP_DATA_MT_HSSPHN_ONLY))
! Dump in HDF only if explicitly requested
                  inquire(file=".DUMP_HS", exist=dump_hs)
                  if (dump_hs) then
                    call write_hdf_hs(a, b, nk, nv(1))
                  endif
#endif
               ENDIF
c
c--->         update with vacuum terms
c
               IF (film .AND. .NOT.odi%d1) THEN
                  CALL cpu_time(time1)
                  CALL hsvac(
     >                 npotmatfile,nmzxyd,n2d,k1d,k2d,k3d,n3d,nmzd,jspd,
     >                 nwdd,matsize,nvd,nv2d,nvac,nmzxy,nmz,nlotot,
     >                 jsp,jspins,nw,vzxy(1,1,1,jsp),vz,evac,bbmat,bmat,
     >                 bkpt,k1,k2,k3,kv2,nv,nq2,ig2,ig,rgphs,zrfs,
     >                 invs2,omtil,area,delz,z1,l_noco,l_J,qss,
     >                 n_size,n_rank,
     X                 a,b,
     <                 nv2)
                  CALL cpu_time(time2)
                  cp_hsvac = cp_hsvac + time2 - time1
               ELSEIF (odi%d1) THEN
c-odim
                  CALL cpu_time(time1)
                  CALL od_hsvac(
     >                 npotmatfile,nmzxyd,n2d,k1d,k2d,k3d,n3d,nmzd,jspd,
     >                 nwdd,matsize,nvd,nv2d,nmzxy,nmz,odi%ig,nlotot,
     >                 jsp,jspins,nw,vzxy(1,1,1,jsp),vz,
     >                 evac,bbmat,bmat,
     >                 bkpt,k1,k2,k3,nv,ig,ig2,sk2,phi2,
     >                 omtil,area,delz,z1,odi%M,odi%mb,odi%m_cyl,
     >                 odi%n2d,nq2,
     >                 n_size,n_rank,zrfs,l_noco,l_J,qss,
     <                 a,b,
     <                 nv2)

                  CALL cpu_time(time2)
                  cp_hsvac = cp_hsvac + time2 - time1

c+odim
               END IF

               IF (gwc==2) THEN
                 CALL gw_eig(
     >                      lmaxd,ntypd,neigd,neigd,natd,nop,nvd,jspd,
     >                 lmd,nbasfcn,llod,nlod,nlotot,invtab,kveclo,ntype,
     >                      mrot,taual,neq,lmax,lnonsph,rmt,omtil,bmat,
     >                      bbmat,us(0,1,jsp),dus(0,1,jsp),uds(0,1,jsp),
     >                      duds(0,1,jsp),ddn(0,1,jsp),invsatnr,invsat,
     >                      ulos(1,1,jsp),uulon(1,1,jsp),dulon(1,1,jsp),
     >                      dulos(1,1,jsp),llo,nlo,l_dulo,lapw_l,nrec,
     >                      l_noco,l_ss,jsp,alph,beta,qss,matsize,a,
     >                      ngopr,jmtd,ntypsy,ntypsd,nlh,nlhd,llh,nmem,
     >                      memd,mlh,clnu,rgphs,ig,nq3,k1d,k2d,k3d,n3d,
     >                      jri,dx,
     >                      vrtot(1,0,1,jsp),vpwtot(1,jsp),
     >                      vr   (1,0,1,jsp),vpw   (1,jsp),
     >                      lda_u(:)%l,vs_mmp(-lmaxb,-lmaxb,1,jsp),
     >                      lmaxb,n_u,odi,ods)
                 GOTO 10
               ENDIF

               IF (l_noco) CLOSE (npotmatfile)
#ifdef DUMP_DATA_OLD
! final H and S, just before the eigenproblem solution
               write(misc_name, '(A, ".bin")') 'Hmat_FULL'
               call write_mat(a, misc_name)
               write(misc_name, '(A, ".bin")') 'Smat_FULL'
               call write_mat(b, misc_name)
#endif               
               CALL cpu_time(time1)

               IF (n_size.NE.1) THEN
                 ndim = CEILING(real(neigd)/n_size)
                 ALLOCATE ( z(nmat,ndim), stat = err )
               ELSE
                 ndim = neigd
                 ALLOCATE ( z(nbasfcn,ndim), stat = err )
               ENDIF
               IF (err.NE.0) THEN
                 WRITE (*,*) 'eigen: error during allocation of the'
                 WRITE (*,*) 'eigenvecs',err,'  size: ',nbasfcn*ndim
                 STOP
               ENDIF
               IF (.not.l_wu) THEN
                 IF (n_size.NE.1) THEN

#if ( defined(CPP_T90) || defined(CPP_EVP) ) 
c
c --->           Parallel diagonalizer
c
                   nblw = nv(jsp)
                   ne = ne_sv
                   CALL chani(
     >                        nmat,n,ndim,neigd,
     >                        n_rank,n_size,SUB_COMM,
     >                        l_zref,nred,
     X                        a,b,k1,k2,k3,
     <                        z,eig,ne)
                   evac_sv(:) = evac(1,:,nw)

#ifndef CPP_HDF
                   IF (n_rank+1.LT.n_size) evac(1,jsp,nw)=999.9
                   IF (l_noco.AND.(n_rank+1.LT.n_size)) 
     +                                  evac(1,jspd,nw) = 999.9
#endif
                   nmat = nmat - nblw + nv(jsp)
                   nblw = -9999
#else 
                   STOP 'eigen: n_size =/= 1  &  no CPP_EVP?'
#endif
                 ELSEIF (gw == 0) THEN
#ifndef CPP_F90

                   CALL zsymsecloc(
     >                       film,nvd,k1,k2,k3,bkpt,nlotot,kveclo,
     >                       ntypd,ntype,nlod,nlo,neq,llo,natd,taual,
     >                       ngopr,nop,invtab,mrot,bmat,nv(jsp),
     >                       nbasfcn,matsize,ndim,
     >                       nmat,bound_lo(nw),bound_up(nw),
     >                       l_zref,l_J,matind,nred,
     >                       a,b,invsat,invsatnr,
     <                       z,eig,ne)
#else
                   CALL franza(
     >                         nbasfcn,matsize,ndim,
     >                         nmat,bound_lo(nw),bound_up(nw),
     >                         (l_zref.and.(nlotot.eq.0)),
     >                         l_J,matind,nred,
     >                         a,b,
     <                         z,eig,ne)
#endif
                   nblw = -9999

                 ENDIF
                 CALL cpu_time(time2)
                 cp_seclr4 = cp_seclr4 + time2 - time1
c
               ELSE
                 CALL aline(
     >                      nvd,lmaxd,ntypd,neigd,natd,nop,jspd,nwdd,
     >                      nlod,matsize,lmd,lmplmd,ntype,neq,lmax,
     >                      rmt,taual,bmat,mrot,ngopr,omtil,film,
     >                      invsat,nlotot,nrec,jsp,nw,el,lnonsph,
     >                      us(0,1,jsp),dus(0,1,jsp),uds(0,1,jsp),
     >                      duds(0,1,jsp),ddn(0,1,jsp),a,b,k1(1,jsp),
     >                      k2(1,jsp),k3(1,jsp),rk(1,jsp),nbasfcn,
     >                      tuu(0,1,1),tud(0,1,1),tdu(0,1,1),
     >                      tdd(0,1,1),ind(0,0,1,1),invtab,invsatnr,
     >                      ulos,uulon,dulon,dulos,llod,llo,nlo,
     >                      l_dulo,lapw_l,l_noco,l_ss,alph,beta,qss,
     >                      odi,ods,bbmat,tuulo(0,-llod,1,1),
     >                      tdulo(0,-llod,1,1),tuloulo(-llod,-llod,1,1),
     >                      mlotot,mlolotot,ello,cp_sph,cp_nsph,
     <                      cp_seclr4,cp_hssphn,bkpt,nv,z,eig,ne)
                   nblw = -9999
               ENDIF
c--->          SECOND VARIATION STEP
               IF (secvar) THEN
c--->            output results from first diagonalization
                  CALL cpu_time(time1)
                  CALL outeig_muff(
     >                             neigd,bound_lo(nw),lpr,film,eonly,
     >                             bkpt,nv(jsp),eig,ne,nblw,nv2(jsp))
                  CALL cpu_time(time2)
                  cp_outeigm = cp_outeigm + time2 - time1
c--->            compute and diagonalize the second variation 
c--->            hamiltonian
                  CALL cpu_time(time1)
                  CALL aline_muff(
     >                            nvd,lmaxd,ntypd,neigd,natd,nop,
     >                            lmd,lmplmd,ntype,neq,lmax,rmt,taual,
     >                            bmat,mrot,ngopr,omtil,invsat,
     >                            nb0,jsp,ne,nbasfcn,
     >                            us,dus,uds,duds,bkpt,rk(1,jsp),
     >                            k1(1,jsp),k2(1,jsp),k3(1,jsp),nv(jsp),
     X                            z,eig)
                  CALL cpu_time(time2)
                  cp_aline = cp_aline + time2 - time1
               END IF
c
c--->         output results
c
               CALL cpu_time(time1)
#ifndef CPP_HDF
               CALL outeig(
     >                     nrec,nvd,neigd,lmaxd,ntypd,jspd,irank,isize,
     >                     bound_lo(nw),wtkpt(nk,nw),lpr,film,
     >                     l_noco,l_ss,l_J,form66,eonly,jsp,nw,
     >                     el(0,1,1,nw),evac(1,1,nw),bkpt,nv,k1,k2,k3,
     >                     z,eig,ne,nblw,nv2,nbasfcn,nmat,ello,nlod,
     >                     nlotot,kveclo,nkpt_l,nkpt(nw),
     <                     eig_l)
#else
#if ( !defined( CPP_INVERSION) )
               IF (.not.l_J) THEN
                  z(:nmat,:ne) = conjg(z(:nmat,:ne))
               ELSE
                  z(:nmat,:ne) = cmplx(0.0,0.0)
               ENDIF
#endif
               eig_l(1,1) = 0.0 ! not used, if hdf is linked
               CALL writeeig(
     >                       nk,jsp,ne,nv(jsp),nmat,k1(:nv(jsp),jsp),
     +                       k2(:nv(jsp),jsp),k3(:nv(jsp),jsp),bkpt,
     >                       wtkpt(nk,nw),eig(:ne),el(0:,:,jsp,nw),
     >                       ello(:,:,jsp),evac(:,jsp,nw),z(:nmat,:ne),
     >                       nlotot,kveclo,n_size,n_rank,neigd)
               IF (l_noco) THEN
                  CALL writebas(
     >                          nk,jspd,nv(jspd),k1(:nv(jspd),jspd),
     >                          k2(:nv(jspd),jspd),k3(:nv(jspd),jspd),
     >                          nlotot,kveclo,el(0:,:,jspd,nw),
     >                          ello(:,:,jspd),evac(:,jspd,nw) )
               ENDIF
#endif
               ! WRITE (*,*) 'pe: ',irank,' wrote ',nrec
               CALL cpu_time(time2)
!#ifdef CPP_MPI
#if ( defined(CPP_T90)  || defined(CPP_HDF) ) 
               IF (n_size.NE.1) evac(1,:,nw) = evac_sv(:)
               !IF (l_noco.AND.(n_size.NE.1)) evac(1,2,nw) = evac_sv
#endif
               cp_outeig = cp_outeig + time2 - time1
               DEALLOCATE ( z )
#ifdef DUMP_SINGLE_KPT
! stop after the 1st k point
            stop 'stopping after 1st k point'
#endif
c
   10       CONTINUE                        ! k-points loop ends
            nb0 = nb0 + ntype
            DEALLOCATE (tuu,tud,tdu,tdd,ind,tuulo,tdulo,tuloulo)
   20    CONTINUE                           ! spin loop ends
!         IF ((n_u.GT.0).AND.(nw.EQ.nwd)) 
         DEALLOCATE( vs_mmp )
         nk0 = nk0 + nkpt(nw)*n_size
         DEALLOCATE (a,b,matind)
   30 CONTINUE                              ! windows loop ends
  333  IF ((gw.eq.2).AND.(gwc==2)) DEALLOCATE ( vpwtot, vrtot )

      IF ( gw.eq.2.AND.(gwc==1) )  THEN        ! go for another round
        OPEN (nu,file='potcoul',form='unformatted',status='old')
!
!       Generate input file abcoeff for subsequent GW calculation
!       28.10.2003 Arno Schindlmayr
!
        WRITE(6,'(A)') 'Info: Write out vxc for GW and vxc.old.'
        WRITE(6,'(A)') 'Info: Write out abcoeff for GW.'
        WRITE(6,'(A)') 'Info: Write out radfun for gw_vxc and GW.'
        OPEN (12,file='vxc.old',form='formatted',status='unknown') ! contains vxc from gw_eig
        OPEN (13,file='vxc',form='formatted',status='unknown')     ! contains vxc from gw_vxc
        OPEN (15,file='abcoeff',form='unformatted',status='unknown',
     &    action='write')
        OPEN (14,file='radfun',form='unformatted',status='unknown')
        WRITE(14) jri(1:ntype)
        gwc=2
        GOTO 999
      ELSE IF ( gw.eq.2.AND.(gwc==2) )  THEN 
        CLOSE (12)
        CLOSE (13)
        CLOSE (14)
        CLOSE (15)
        STOP 'GW finished'
      ENDIF
      DEALLOCATE ( vpw,vzxy,vz,vr,vr0 )

      IF (irank.EQ.0) THEN
        CALL outtime('   tlmplm:',cp_tlmplm)
        CALL outtime('   generation of planewaves (apws):',cp_apws)
        CALL outtime('   hsint:',cp_hsint)
        WRITE (2,*)'  hssph(n):'
        IF (secvar) THEN
           CALL outtime('      hssph (total):',cp_hssphn)
           CALL outtime('      spherical:',cp_sph)
        ELSE
           CALL outtime('      spherical:',cp_sph)
           CALL outtime('      non-spherical:',cp_nsph)
           CALL outtime('      hssphn (total):',cp_hssphn)
        END IF
        IF (film) THEN
           CALL outtime('   hsvac:',cp_hsvac)
        END IF
        CALL outtime('   seclr4:',cp_seclr4)
        IF (secvar) THEN
          CALL outtime('   outeig_muff:',cp_outeigm)
          CALL outtime('   second variation step (aline_muff)',cp_aline)
        END IF
        CALL outtime('   outeig:',cp_outeig)
      ENDIF

#if (defined(CPP_MPI) && !defined(CPP_T90))
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
#ifdef CPP_HDF
      CALL closeeig()
#else
      IF (l_J) THEN
        IF (irank==0) CLOSE (66)
      ELSE
        CLOSE (66)
      ENDIF
#endif

      IF (form66) STOP 'form66.complete'
      IF (l_f.AND.(irank.EQ.0)) THEN
        CLOSE (28)
        IF (jspins.EQ.2) CLOSE (38)
      ENDIF

      DEALLOCATE ( nv,nv2,k1,k2,k3,bkpt,ddn,uloulopn,eig,rk,ello )
      DEALLOCATE ( el,bound_lo,bound_up,us,uds,dus,duds,ulos,dulos )
      DEALLOCATE ( evac,uulon,dulon,kveclo )

      END SUBROUTINE eigen
      END MODULE m_eigen
