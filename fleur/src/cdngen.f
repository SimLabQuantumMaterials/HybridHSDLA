      MODULE m_cdngen
      CONTAINS
      SUBROUTINE cdngen(
     > nrhomfile,npotmatfile,noinpfile,kcrel,
     > cdinf,dos,frcor,slice,ctail,ndir,vacdos,integ,layers,izlay,
     > kk,e1s,e2s,nnne,pallst,layerd,neigd,nkptd,nv2d,nbasfcn,
     > ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd,
     > nvd,nmzd,k1d,k2d,k3d,kq1d,kq2d,kq3d,msh,
     > nop,n2d,natd,nwdd,nmzxyd,nspd,nn3d,nstd,irecl0,
     > namat,igfft,pgfft,ig,sk3,nstr,nstr2,kv2,enmix,
     > amat,bmat,mx1,mx2,mx3,film,symor,zrfs,invs2,
     > invs,ngopr,volint,vol,volmts,zatom,jri,l_soc,lda_u,n_u,
     > jspins,nwd,nkpt,ntype,nq3,nmz,nq2,nvac,nmzxy,
     > nsymt,clnu,llh,mlh,nlh,nmem,ntypsy,pos,taual,
     > rmt,rmsh,dx,neq,lmax,z1,omtil,area,delz,theta,phi,
     > kq1_fft,kq2_fft,kq3_fft,nq3_fft,kmxq2_fft,kmxq_fft,kv3,
     > nop2,mrot,zelec,tau,ig2,rgphs,irank,isize,ncst,
     > lchange,lchg_v,w_iks,l_f,l_geo,bbmat,invarop,invarind,
     > multab,invtab,invsat,invsatnr,d_wgn,nlod,llod,nlo,llo,
     > l_dulo,ulo_der,lapw_l,llochg,skiplo,nlol,lo1l,gw,
     > lepr,sigma,tworkf,nstars,nstm,starcoeff,locx,locy,
     > fname,sk2,phi2,odi,ods,
     X l_noco,l_ss,qss,sso_opt,l_mperp,l_relax,l_constr,mix_b,
     X ello0,el0,evac0,force,alph,beta,b_con,
     < seigc,e1_dos,e2_dos,efermi,sig_dos)
c
c     *****************************************************
c     flapw7 charge density generator
c                                c.l.fu
c     correction for nwd>1: semic must be set false for nw=nwd
c     furthermore ch.d. arrays initialized to zero INSIDE nw
c     loop; only last window stored; shift do 10
c     aug 90; r.p.
c     modifications for slicing by r.p. and h.e. (1/91)
c     error in reading 66 for spin-polarized case:
c     read(66) nkpt was read for jspin>1 also
c     Corrected mai 93 r.p.
c     *****************************************************
C
      USE m_constants, ONLY : pimach
      USE m_umix
      USE m_prpqfftmap
      USE m_cdnval
      USE m_loddop
      USE m_wrtdop
      USE m_cdntot
      USE m_cdnovlp
      USE m_qfix
      USE m_cputime
      USE m_outtime
      USE m_rwnoco
      USE m_od_types, ONLY : od_inp, od_sym

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nrhomfile,npotmatfile,noinpfile,irecl0
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd
      INTEGER, INTENT (IN) :: nvd,nmzd,k1d,k2d,k3d,kq1d,kq2d,kq3d,msh
      INTEGER, INTENT (IN) :: nop,n2d,natd,nwdd,nmzxyd,nspd,nn3d,nstd
      INTEGER, INTENT (IN) :: layerd,neigd,nkptd,nwd,nv2d,nbasfcn
      INTEGER, INTENT (IN) :: jspins,ntype,nq3,nmz,nq2,nvac,nmzxy,nsymt
      INTEGER, INTENT (IN) :: layers,kk,nnne,mx1,mx2,mx3,irank,isize
      INTEGER, INTENT (IN) :: llod,nlod,lepr,n_u,gw
      LOGICAL, INTENT (IN) :: ctail,dos,frcor,slice,cdinf
      LOGICAL, INTENT (IN) :: vacdos,integ,pallst,l_f,l_soc
      LOGICAL, INTENT (INOUT) :: l_noco,l_ss,l_mperp,l_constr
      LOGICAL, INTENT (IN) :: film,symor,zrfs,invs2,invs
      REAL,    INTENT (IN) :: e1s,e2s,theta,phi,sigma
      REAL,    INTENT (IN) :: delz,area,omtil,z1,volint,vol
      REAL,    INTENT (IN) :: e1_dos,e2_dos,efermi,sig_dos
      REAL,    INTENT (INOUT) :: mix_b
      REAL,    INTENT (OUT) :: seigc
      INTEGER kcrel,ndir,nop2
c+dw
      REAL,    INTENT (IN) :: tworkf,locx(2),locy(2)
      INTEGER, INTENT (IN) :: nstars,nstm
      LOGICAL, INTENT (IN) :: starcoeff
c-dw
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: d_wgn(-3:3,-3:3,3,nop)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),lmax(ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: nlh(ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),ngopr(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (IN) :: multab(nop,nop),invtab(nop)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd),nkpt(nwdd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: nstr2(n2d),nstr(n3d),kv2(2,n2d),kv3(3,n3d)
      INTEGER, INTENT (IN) :: igfft(0:nn3d-1,2),izlay(layerd,2)
      INTEGER, INTENT (IN) :: kq1_fft(nwdd),kq2_fft(nwdd),kq3_fft(nwdd)
      INTEGER, INTENT (IN) :: nq3_fft(nwdd),kmxq_fft(nwdd),ncst(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),lapw_l(ntypd)
      INTEGER, INTENT (IN) :: skiplo(ntypd,jspd),lda_u(ntypd)
      INTEGER, INTENT (IN) :: nlol(0:llod,ntypd),lo1l(0:llod,ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL, INTENT    (IN) :: sk3(n3d),tau(3,nop),volmts(ntypd)
      REAL, INTENT    (IN) :: amat(3,3),bmat(3,3),bbmat(3,3)
      REAL, INTENT    (IN) :: pos(3,natd),dx(ntypd),zatom(ntypd)
      REAL, INTENT    (IN) :: rmsh(jmtd,ntypd)
      REAL, INTENT    (IN) :: rmt(ntypd),taual(3,natd)
      REAL, INTENT    (IN) :: pgfft(0:nn3d-1)
      REAL, INTENT    (IN) :: w_iks(neigd,nkptd,jspd),zelec(nwdd)
      REAL, INTENT    (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL, INTENT (INOUT) :: alph(ntypd),beta(ntypd),b_con(2,ntypd)
      REAL, INTENT (INOUT) :: el0(0:lmaxd,ntypd,jspd,nwdd),qss(3)
      REAL, INTENT (INOUT) :: evac0(2,jspd,nwdd)
      REAL, INTENT (INOUT) :: ello0(nlod,ntypd,jspd)
      REAL, INTENT (INOUT) :: force(3,ntypd,jspd)
      REAL, INTENT    (IN) :: enmix(jspd,nwdd)
      LOGICAL, INTENT (IN) :: lchange(0:lmaxd,ntypd,jspd,nwdd)
      LOGICAL, INTENT (IN) :: lchg_v(2,jspd,nwdd),l_geo(ntypd)
      LOGICAL, INTENT (IN) :: llochg(nlod,ntypd,jspd),l_dulo(nlod,ntypd)
      LOGICAL, INTENT (INOUT) :: l_relax(ntypd),sso_opt(2) 
      CHARACTER*2, INTENT(IN) :: namat(0:103)
      CHARACTER*12,INTENT(IN) :: fname(3)
c-odim
      REAL,    INTENT (IN) :: sk2(n2d),phi2(n2d)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL fix,qtot,scor,seig,smom,stot,sval,time1,time2,dummy,pi,sfp
      REAL slmom,slxmom,slymom,sum,thetai,phii
      INTEGER iter,ivac,j,jspin,jspmax,k,n,nt,nw,nrec,ieig,ikpt
      INTEGER lmd,llpd,ityp,nk0,ilayer,urec,itype,iatom
      LOGICAL semic,l_relax_any,exst,n_exist
C     ..
C     .. Local Arrays ..
      REAL stdn(ntypd,jspd),svdn(ntypd,jspd),alpha_l(ntypd),
     +     rh(msh,ntypd,jspd),vr0(jmtd,ntypd,jspd),qint(ntypd,jspd)
      REAL chmom(ntypd,jspd),clmom(3,ntypd,jspd),cp_time(9)
      CHARACTER*8 dop,iop,name(10)
      INTEGER,ALLOCATABLE :: igq_fft(:)
      REAL   ,ALLOCATABLE :: vz(:,:,:),vr(:,:,:,:)
      REAL   ,ALLOCATABLE :: rht(:,:,:),rho(:,:,:,:)
      COMPLEX,ALLOCATABLE :: vpw(:,:),vzxy(:,:,:,:)
      COMPLEX,ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:)
c---> pk non-collinear
      REAL    rhoint,momint,alphdiff(ntypd)
      INTEGER kmxq2_fft(nwdd),igq2_fft(0:kq1d*kq2d-1)
      REAL   ,ALLOCATABLE :: qvac(:,:,:,:),qvlay(:,:,:,:,:),qis(:,:,:)
      COMPLEX,ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:),qa21(:)
c---> pk non-collinear
c--- J<
      REAL thetaJ
      REAL M(ntypd)
      INTEGER qnumb,nmagn,nsh,mtypes
      INTEGER magtype(ntypd),nmagtype(ntypd)
      LOGICAL l_J,l_disp
      LOGICAL l_magn(ntypd)
c--- J>
C     ..
C     .. External Subroutines ..
      EXTERNAL cored,coredr,m_perp
C     ..
      pi  = pimach()
      sfp = 2* sqrt( pimach() )
      lmd  = lmaxd*(lmaxd+2)
      llpd = (lmaxd*(lmaxd+3))/2

c---> timing
      IF (irank.EQ.0) WRITE (2,8005)
 8005 FORMAT ('CHARGE DENSITY PART (cdngen):')
      DO j=1,9
        cp_time(j) = 0.0       ! 1=cdnval ; 2=pwden ; 3=vacden ; 4=abcof  
      ENDDO                    ! 5=cdinf  ; 6=rhomt ; 7=rhonmt ; 8=mtsum ; 9 = communication
c
c Read Potential and keep only vr(:,0,:,:) and vz
c
      ALLOCATE(vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd),
     +       vz(nmzd,2,jspd),vr(jmtd,0:nlhd,ntypd,jspd))
      OPEN (8,file='pottot',form='unformatted',status='old')
      REWIND (8)
      CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,8,natd,neq,
     <            iop,dop,iter,vr,vpw,vz,vzxy,name)
      CLOSE (8)
      DO jspin = 1,jspins 
        DO n = 1, ntype
          DO j = 1,jri(n)
            vr0(j,n,jspin) = vr(j,0,n,jspin) 
          ENDDO
        ENDDO 
      ENDDO 
      DEALLOCATE ( vpw,vzxy )
      ALLOCATE ( qpw(n3d,jspd),rhtxy(nmzxyd,odi%n2d-1,2,jspd) )
      ALLOCATE ( rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd) )
c
c Read in input density
c
      nt = 71
      IF ( (.NOT. l_noco) .AND. irank.EQ.0) THEN
        OPEN (nt,file='cdn1',form='unformatted',status='old')
        REWIND (nt)
        CALL loddop(
     >              jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >              jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >              nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <              iop,dop,iter,rho,qpw,rht,rhtxy,name)
      ENDIF
c
c     *** exchanged the order of loops (first windows than spin)   ***
c     *** in order to read the data from file eig (unit 66) in the ***
c     *** same order as it is written in eigen. Otherwise the      ***
c     *** program doesn't work for spinpolarised calculations with ***
c     *** more than one window.     p.kurz ***
c
c====> loop over windows
c
      IF (irank.EQ.0) 
     +         OPEN (40,file='enpara',form='formatted',status='unknown')
      nk0 = 0
      nrec = 0
      ALLOCATE ( cdom(n3d),cdomvz(nmzd,2),cdomvxy(nmzxyd,odi%n2d-1,2) )
      ALLOCATE ( qvlay(neigd,layerd,2,nkptd,jspd),qa21(ntypd) )
      ALLOCATE ( qvac(neigd,2,nkptd,jspd),qis(neigd,nkptd,jspd) )
      ALLOCATE ( igq_fft(0:kq1d*kq2d*kq3d-1) )
c
      DO 110 nw = 1,nwd
!
!--->    initialize density arrays with zero
!
         IF ((nw.EQ.1).OR.(irank.NE.0)) THEN
           qa21(:) = cmplx(0.0,0.0)
           rho(:,:,:,:) = 0.0
           qpw(:,:) = cmplx(0.0,0.0)
           cdom(:) =  cmplx(0.0,0.0) 
           IF (film) THEN
              rht(:,:,:) = 0.0
              cdomvz(:,:) = cmplx(0.0,0.0)
              rhtxy(:,:,:,:) = cmplx(0.0,0.0)
              cdomvxy(:,:,:) = cmplx(0.0,0.0)
            END IF
         ENDIF
!
!--->    Set up pointer for backtransformation of from g-vector in
!        positive domain fof carge density fftibox into stars
!        In principle this can also be done in main program once.
!        It is done here to save memory.
!
         CALL prp_qfft_map(
     >                     n3d,kq1d,kq2d,kq3d,nop,
     >                     kv3,nop2,mrot,film,
     >                     kq1_fft(nw),kq2_fft(nw),kq3_fft(nw),
     >                     nq3_fft(nw),kmxq2_fft(nw),kmxq_fft(nw),
     <                     igq2_fft,igq_fft)

         semic = .false.
         IF (nw.NE.nwd) semic = .true.
         IF (film) THEN
            qvac(:,:,:,:) = 0.0
            qvlay(:,:,:,:,:) = cmplx(0.0,0.0)
         ENDIF
!
!--->    LDA+U: initialise density-matrix if needed
!
         IF (n_u.GT.0) THEN
            ALLOCATE (  n_mmp(-3:3,-3:3,n_u,jspins) )
            n_mmp(:,:,:,:) = cmplx(0.0,0.0)
         ELSE
            ALLOCATE ( n_mmp(-3:-3,-3:-3,1,jspins))
         ENDIF

!
!--->    in a non-collinear calcuation where the off-diagonal part of
!        density matrix in the muffin-tins is calculated, the a- and 
!        b-coef. for both spins are needed at once. Thus, cdnval is only
!        called once and both spin directions are calculated in a single
!        go.
!
         jspmax = jspins
         IF (l_mperp) jspmax = 1
         DO 100 jspin = 1,jspmax 
           CALL cpu_time(time1)
           IF (l_noco) nrec = nk0
c+fo
           IF (l_f.AND.(jspin.EQ.1)) THEN
             INQUIRE (file=fname(2),exist=exst)
             IF (.NOT.exst) STOP 'cdnval & force_a21: file tmat missing'
             OPEN (28,file=fname(2),access='direct',recl=irecl0)
           ENDIF
           IF (l_f.AND.(jspin.EQ.2)) THEN
             INQUIRE (file=fname(3),exist=exst)
             IF (.NOT.exst) STOP 'cdnval & force_a21: file tmas missing'
             OPEN (38,file=fname(3),access='direct',recl=irecl0)
           END IF
c-fo
        CALL cdnval(
     >   66,npotmatfile,nkpt(nw),jspin,nw,semic,slice,l_noco,
     >   l_ss,l_mperp,cdinf,dos,ndir,vol,volint,volmts,enmix,sk3,vacdos,
     >   integ,layers,izlay,kk,e1s,e2s,nnne,pallst,layerd,neigd,nkptd,
     >   nv2d,ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd,nvd,nmzd,k1d,
     >   k2d,k3d,kq1d,kq2d,kq3d,irank,isize,nbasfcn,nop,n2d,natd,nwdd,
     >   nmzxyd,nspd,lmd,llpd,
     >   kq1_fft(nw),kq2_fft(nw),kq3_fft(nw),nq3_fft(nw),kmxq_fft(nw),
     >   igq_fft,igfft,pgfft,zelec,w_iks,vr0,
     >   vr(1,0,1,jspin),vz(1,1,jspin),lepr,e1_dos,e2_dos,efermi,
     >   sig_dos,lchange,lchg_v,l_f,l_geo,bbmat,sfp,invarop,invarind,
     >   multab,invtab,invsat,invsatnr,ntype,film,zrfs,invs,
     >   symor,nvac,nmz,nmzxy,l_soc,jspins,delz,area,omtil,z1,amat,bmat,
     >   taual,pos,nop2,llh,nmem,mlh,clnu,nsymt,ngopr,nlh,ntypsy,mrot,
     >   tau,nstr,nstr2,ig,ig2,rgphs,mx1,mx2,mx3,kv2,kv3,nq3,nq2,rmt,dx,
     >   rmsh,jri,lmax,neq,nlod,llod,nlo,llo,llochg(1,1,jspin),
     >   lda_u,n_u,skiplo(1,jspin),nlol,lo1l,l_dulo,ulo_der,lapw_l,
     >   tworkf,nstars,nstm,starcoeff,locx,locy,msh,nstd,
     >   odi,ods,sk2,phi2,ncst,zatom,d_wgn,alph,beta,qss,
     X   nrec,n_mmp(-3,-3,1,jspin),cp_time,ello0,el0,evac0,force,
     X   qpw,rhtxy,rho,rht,qvac,qvlay,cdom,cdomvz,cdomvxy,qa21,
     <   chmom,clmom,qis)
           CALL cpu_time(time2)
           cp_time(1) = cp_time(1) + time2 - time1 ! 1=cdnval
c+fo
           IF (l_f.AND.(jspin.EQ.1)) CLOSE (28)
           IF (l_f.AND.(jspin.EQ.2)) CLOSE (38)
c-fo
 100     CONTINUE
c-lda+U
         IF ( (n_u.GT.0).AND.(nw.EQ.nwd) ) THEN
           IF (irank.EQ.0) CALL u_mix(
     >                                n_u,jspins,n_mmp)
         ENDIF
         DEALLOCATE (  n_mmp )
c-lda-U
         nk0 = nk0 + nkpt(nw)
 110  CONTINUE
      DEALLOCATE (vz,vr)
c+t3e      
      IF (irank.EQ.0) THEN
c-t3e
      CLOSE (40)
      WRITE (2,*)'  cdnval:'
      CALL outtime('      pwden:',cp_time(2))
      IF (film) CALL outtime('      vacden:',cp_time(3))
      CALL outtime('      readin:',cp_time(3))
      CALL outtime('      abcof:',cp_time(4))
      IF (dos .OR. cdinf) THEN
         CALL outtime('      cdninf:',cp_time(5))
      ENDIF
      CALL outtime('      rhomt:',cp_time(6))
      CALL outtime('      rhonmt:',cp_time(7))
      CALL outtime('      sum of charge in spheres:',cp_time(8))
      IF (isize > 0) CALL outtime('      communication :',cp_time(9))
      CALL outtime('      cdnval (total):',cp_time(1))

      CALL cdntot(
     >            k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,
     >            nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,
     >            ntype,neq,volmts,taual,z1,vol,volint,
     >            symor,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >            nstr,kv3,delz,jri,dx,rmsh,invtab,
     >            qpw,rho,rht,odi,
     <            qtot,dummy)
c
c---> changes
c
      ENDIF ! irank = 0
      IF (kcrel.EQ.0) THEN
         seigc = 0.
!
! Generate input file ecore for subsequent GW calculation
! 11.2.2004 Arno Schindlmayr
!
         IF ((gw.eq.1).AND.(irank.EQ.0)) THEN
            OPEN (15,file='ecore',status='unknown',
     &                     action='write',form='unformatted')
         ENDIF

         DO 131 jspin = 1,jspins
            IF ((jspins.EQ.2).AND.(irank.EQ.0)) THEN
               DO 111 n = 1,ntype
                  svdn(n,jspin) = rho(1,0,n,jspin)/
     +                            (sfp*rmsh(1,n)*rmsh(1,n))
 111           CONTINUE
            END IF
c
c     block 1 unnecessary for slicing: begin
            IF (.NOT.slice) THEN
c     ---> add in core density
              IF (irank.EQ.0) THEN
               CALL cored(
     >                 jspins,jspin,ntype,neq,zatom,rmsh,dx,jri,film,
     >                 frcor,ncst,rmt,rho,jmtd,jspd,msh,
     >                 nlhd,nstd,ntypd,
     >                 vr0(1,1,jspin),gw,
     <                 qint,rh,seig)
               seigc = seigc + seig
               IF (jspins.EQ.2) THEN
                  DO 121 n = 1,ntype
                     stdn(n,jspin) = rho(1,0,n,jspin)/
     +                               (sfp*rmsh(1,n)*rmsh(1,n))
 121              CONTINUE
               ENDIF
              ENDIF  ! irank = 0
c     ---> add core tail charge to qpw
              IF ((l_noco).AND.(irank.EQ.0)) THEN
c---> pk non-collinear
c--->           add the coretail-charge to the constant interstitial
c--->           charge (star 0), taking into account the direction of
c--->           magnetisation of this atom
                IF (jspin .EQ. 2) THEN
                  DO ityp = 1,ntype
                     rhoint  = (qint(ityp,1) + qint(ityp,2))
     /                          /volint/jspins/2.0
                     momint  = (qint(ityp,1) - qint(ityp,2))
     /                          /volint/jspins/2.0
c--->                rho_11
                     qpw(1,1) = qpw(1,1)
     +                  + rhoint + momint*cos(beta(ityp))
c--->                rho_22
                     qpw(1,2) = qpw(1,2)
     +                  + rhoint - momint*cos(beta(ityp))
c--->                real part rho_21
                     cdom(1) = cdom(1) + cmplx(0.5*momint
     *                         *cos(alph(ityp))*sin(beta(ityp)),0.0)
c--->                imaginary part rho_21
                     cdom(1) = cdom(1) + cmplx(0.0,-0.5*momint
     *                         *sin(alph(ityp))*sin(beta(ityp)))
                  ENDDO
                ENDIF
c---> pk non-collinear
              ELSEIF (ctail) THEN
                CALL cdnovlp(irank,isize,
     >                   memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop,
     >                   n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh,
     >                   ntype,mx3,nq2,nq3,ntypsy,lmax,nmzxy,
     >                   nmz,nvac,neq,bmat,taual,kv3,sk3,ig,kv2,nstr2,
     >                   nstr,film,zrfs,z1,delz,omtil,rmt,dx,rmsh,
     >                   jri,mlh,clnu,llh,nlh,nmem,mrot,tau,symor,
     >                   jspin,rh,invs,ncst,invtab,
     >                   odi,ods,amat,ig2,sk2,phi2,vol,
     X                   qpw,rhtxy,rho,rht)
              ELSEIF (irank.EQ.0) THEN
                  DO ityp = 1,ntype
                     qpw(1,jspin) = qpw(1,jspin)
     +                            + qint(ityp,jspin)/jspins/volint
                  ENDDO
              END IF
c     block 1 unnecessary for slicing: end
            END IF
c
 131     ENDDO

         IF ((gw.eq.1).AND.(irank.EQ.0)) CLOSE(15)
      ELSE
c relativistic core implementation : kcrel.eq.1
         seigc = 0.
         IF ((jspins.EQ.2).AND.(irank.EQ.0)) THEN
            DO jspin = 1,jspins
               DO n = 1,ntype
                  svdn(n,jspin) = rho(1,0,n,jspin)/
     +                            (sfp*rmsh(1,n)*rmsh(1,n))
               END DO
            END DO
         END IF
c
c     block 1 unnecessary for slicing: begin
         IF (.NOT.slice) THEN
c     ---> add in core density
           IF (irank.EQ.0) THEN
            CALL coredr(jspins,ntype,neq,zatom,rmsh,dx,jri,seig,film,
     >                  ncst,rmt,rho,jmtd,jspd,msh,nlhd,nstd,ntypd,
     >                  vr0,
     <                  qint,rh)
            seigc = seigc + seig
            IF (jspins.EQ.2) THEN
               DO jspin = 1,jspins
                  DO n = 1,ntype
                     stdn(n,jspin) = rho(1,0,n,jspin)/
     +                               (sfp*rmsh(1,n)*rmsh(1,n))
                  END DO
               END DO
            END IF
           ENDIF
            IF ((l_noco).AND.(irank.EQ.0)) THEN
c---> pk non-collinear
c--->          add the coretail-charge to the constant interstitial
c--->          charge (star 0), taking into account the direction of
c--->          magnetisation of this atom
               DO ityp = 1,ntype
                  rhoint  = (qint(ityp,1) + qint(ityp,2))
     /                          /volint/jspins/2.0
                  momint  = (qint(ityp,1) - qint(ityp,2))
     /                          /volint/jspins/2.0
c--->             rho_11
                  qpw(1,1) = qpw(1,1)
     +                 + rhoint + momint*cos(beta(ityp))
c--->             rho_22
                  qpw(1,2) = qpw(1,2)
     +                  + rhoint - momint*cos(beta(ityp))
c--->             real part rho_21
                  cdom(1) = cdom(1) + cmplx(0.5*momint
     *                         *cos(alph(ityp))*sin(beta(ityp)),0.0)
c--->             imaginary part rho_21
                  cdom(1) = cdom(1) + cmplx(0.0,-0.5*momint
     *                         *sin(alph(ityp))*sin(beta(ityp)))
               ENDDO
c---> pk non-collinear
            ELSE
               DO jspin = 1,jspins
                  IF (ctail) THEN
c+gu hope this works as well
                     CALL cdnovlp(irank,isize,
     >                   memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop,
     >                   n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh,
     >                   ntype,mx3,nq2,nq3,ntypsy,lmax,nmzxy,
     >                   nmz,nvac,neq,bmat,taual,kv3,sk3,ig,kv2,nstr2,
     >                   nstr,film,zrfs,z1,delz,omtil,rmt,dx,rmsh,
     >                   jri,mlh,clnu,llh,nlh,nmem,mrot,tau,symor,
     >                   jspin,rh(1,1,jspin),invs,ncst,invtab,
     >                   odi,ods,amat,ig2,sk2,phi2,vol,
     X                   qpw,rhtxy,rho,rht)
                  ELSEIF (irank.EQ.0) THEN
                     DO ityp = 1,ntype
                        qpw(1,jspin) = qpw(1,jspin)
     +                       + qint(ityp,jspin)/jspins/volint
                     ENDDO
                  END IF
               END DO
            ENDIF
c     block 1 unnecessary for slicing: end
         END IF
c end relativistic core
      END IF
      IF (irank.EQ.0) THEN
c     block 2 unnecessary for slicing: begin
      IF (.NOT.slice) THEN
         CALL qfix(
     >             k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >             nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,nmzxy,n2d,
     >             ntype,neq,volmts,taual,z1,vol,volint,nq2,invtab,
     >             symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >             nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     >             qpw,rhtxy,rho,rht,odi,
     <             fix)
c---> pk non-collinear
         IF (l_noco) THEN
c--->    fix also the off-diagonal part of the density matrix
         DO k = 1,nq3
            cdom(k) = fix*cdom(k)
         ENDDO
         IF (film) THEN
            DO ivac = 1,2
               DO j = 1,nmz
                  cdomvz(j,ivac) = fix*cdomvz(j,ivac)
               ENDDO
               DO k = 2,odi%nq2
                  DO j = 1,nmzxy
                     cdomvxy(j,k-1,ivac) = fix*cdomvxy(j,k-1,ivac)
                  ENDDO
               ENDDO
            ENDDO
         END IF
         ENDIF
c---> pk non-collinear

c     ---> spin densities at the nucleus
c     ---> and magnetic moment in the spheres
         IF (jspins.EQ.2) THEN
            WRITE (6,FMT=8000)
            WRITE (16,FMT=8000)
            DO 150 n = 1,ntype
               sval = svdn(n,1) - svdn(n,jspins)
               stot = stdn(n,1) - stdn(n,jspins)
               scor = stot - sval
               WRITE (6,FMT=8010) n,stot,sval,scor,svdn(n,1),stdn(n,1)
               WRITE (16,FMT=8010) n,stot,sval,scor,svdn(n,1),stdn(n,1)
  150       CONTINUE
            IF (l_mperp) THEN
              ! angles in nocoinp file are (alph-alphdiff)
              iatom= 1
              DO n= 1,ntype
                IF (l_ss) THEN
                  alphdiff(n)= 2.*pi*(  qss(1)*taual(1,iatom)
     +                                + qss(2)*taual(2,iatom)
     +                                + qss(3)*taual(3,iatom) )
                ELSE
                  alphdiff(n)= 0.
                ENDIF
                iatom= iatom + neq(n)
              ENDDO
            ENDIF
            WRITE (6,FMT=8020)
            WRITE (16,FMT=8020)
            DO 160 n = 1,ntype
               smom = chmom(n,1) - chmom(n,jspins)
               WRITE (6,FMT=8030) n,smom, (chmom(n,j),j=1,jspins)
               WRITE (16,FMT=8030) n,smom, (chmom(n,j),j=1,jspins)
               IF (l_mperp) THEN
c--->             calculate the perpendicular part of the local moment
c--->             and relax the angle of the local moment or calculate
c--->             the constraint B-field.
                  CALL m_perp(
     >                 jspd,ntypd,jmtd,n,l_relax,l_constr,mix_b,
     >                 jri,dx,rmsh,vr0,volmts,chmom,qa21,alphdiff,
     X                 alph,beta,b_con)
               ENDIF
  160       CONTINUE

c--->       save the new nocoinp file if the dierctions of the local
c--->       moments are relaxed or a constraint B-field is calculated.
            l_relax_any = .false.
            iatom = 1
            DO itype = 1,ntype
              l_relax_any = l_relax_any.OR.l_relax(itype)
            ENDDO 
            IF (l_relax_any.OR.l_constr) THEN
               IF (.not. l_mperp) THEN
                 STOP 
     &            'cdngen: (l_relax_any.OR.l_constr).AND.(.NOT.l_mperp)' 
               ENDIF 
               DO itype = 1,ntype 
                 IF ( l_ss ) THEN
                   alpha_l(itype) = alph(itype) - alphdiff(itype) 
                   DO WHILE (alpha_l(n) > +pi)
                     alpha_l(n)= alpha_l(n) - 2.*pi
                   ENDDO
                   DO WHILE (alpha_l(n) < -pi)
                     alpha_l(n)= alpha_l(n) + 2.*pi
                   ENDDO
                 ELSE
                   alpha_l(itype) = alph(itype)
                 ENDIF
               ENDDO

               OPEN (noinpfile,file='nocoinp',form='formatted',
     +                                            status='old')
               REWIND (noinpfile)
               CALL rw_noco(
     >                      'W',noinpfile,ntype,l_J,l_soc,
     X                      alpha_l,beta,l_relax,b_con,
     X                      l_ss,l_mperp,l_constr,mix_b,qss,sso_opt,
     <                      thetaJ,l_magn,nmagn,M,mtypes,magtype,
     <                      nmagtype,nsh,l_disp)
               CLOSE (noinpfile)
            ENDIF

            IF (l_soc) THEN
              thetai = theta
              phii   = phi
              WRITE (6,FMT=9020)
              WRITE (16,FMT=9020)
              DO n = 1,ntype
                 IF (l_noco) THEN
                   thetai = beta(n)
                   phii =   alph(n)
                 ENDIF
c
c magn. moment(-)
c
                 slxmom = clmom(1,n,1)+clmom(1,n,2)
                 slymom = clmom(2,n,1)+clmom(2,n,2)
                 slmom =  clmom(3,n,1)+clmom(3,n,2)
c
c  rotation: orbital moment || spin moment (extended to incude phi - hopefully)
c
                 slmom   = cos(thetai)*slmom + sin(thetai)*
     +                     (cos(phii)*slxmom + sin(phii)*slymom)
                 clmom(3,n,1) = cos(thetai)*clmom(3,n,1) + sin(thetai)*
     +                (cos(phii)*clmom(1,n,1) + sin(phii)*clmom(2,n,1))
                 clmom(3,n,2) = cos(thetai)*clmom(3,n,2) + sin(thetai)*
     +                (cos(phii)*clmom(1,n,2) + sin(phii)*clmom(2,n,2))
                 WRITE (6,FMT=8030) n,slmom,(clmom(3,n,j),j=1,2)
                 WRITE (16,FMT=8030) n,slmom,(clmom(3,n,j),j=1,2)
c                WRITE (16,FMT=8030) n,slxmom,(clmom(1,n,j),j=1,2)
c                WRITE (16,FMT=8030) n,slymom,(clmom(2,n,j),j=1,2)
              END DO
            END IF
         END IF
c     block 2 unnecessary for slicing: end
      END IF
 9020 FORMAT (/,/,10x,'orb. magnetic moments in the spheres:',/,10x,
     +       'type',t22,'moment',t33,'spin-up',t43,'spin-down')
 8000 FORMAT (/,/,10x,'spin density at the nucleus:',/,10x,'type',t25,
     +       'total',t42,'valence',t65,'core',t90,
     +       'majority valence and total density',/)
 8010 FORMAT (i13,2x,3e20.8,5x,2e20.8)
 8020 FORMAT (/,/,2x,'-->  magnetic moments in the spheres:',/,2x,
     +       'mm -->   type',t22,'moment',t33,'spin-up',t43,'spin-down')
 8030 FORMAT (2x,'--> mm',i8,2x,3f12.5)
c
      IF (slice) THEN
         OPEN (20,file='cdn_slice',form='unformatted',status='unknown')
         CALL wrtdop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >               jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >               invs,invs2,film,delz,z1,dx,rmt,zatom,
     >               nlh,jri,ntypsd,ntypsy,namat,20,
     >               'slice   ','density ',iter,rho,qpw,rht,rhtxy,name)
         IF (l_noco) THEN
            WRITE (20) (cdom(k),k=1,nq3)
            IF (film) THEN
               WRITE (20) ((cdomvz(j,ivac),j=1,nmz),ivac=1,nvac)
               WRITE (20) (((cdomvxy(j,k-1,ivac),j=1,nmzxy),k=2,odi%nq2)
     +                                              ,ivac=1,nvac)
            ENDIF
         ENDIF
         CLOSE(20) 
         STOP 'slice OK'
      END IF
c
      IF (l_noco) THEN
c---> pk non-collinear
c--->    write output density matrix on file rhomat_out
c--->    first the diagonal elements of the density matrix
         OPEN (nrhomfile,FILE='rhomat_out',FORM='unformatted',
     +         STATUS='unknown')
         CALL wrtdop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >               jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >               invs,invs2,film,delz,z1,dx,rmt,zatom,
     >               nlh,jri,ntypsd,ntypsy,namat,nrhomfile,
     >               'output  ','density ',iter,rho,qpw,rht,rhtxy,name)
c--->    and then the off-diagonal part
         WRITE (nrhomfile) (cdom(k),k=1,nq3)
         IF (film) THEN
            WRITE (nrhomfile) ((cdomvz(j,ivac),j=1,nmz),ivac=1,nvac)
            WRITE (nrhomfile) (((cdomvxy(j,k-1,ivac),j=1,nmzxy),
     +                                    k=2,odi%nq2),ivac=1,nvac)
         ENDIF
         CLOSE (nrhomfile)
c---> pk non-collinear
      ELSE
c      ----> write output density on unit 71
         CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt,
     >            'output  ','density ',iter,rho,qpw,rht,rhtxy,name)
         CLOSE (nt)
      ENDIF
c+t3e
      ENDIF
c+dw added for other pe than pe=0
      DEALLOCATE (qvac,qvlay,qis,cdom,cdomvz,cdomvxy,qa21)
      DEALLOCATE (qpw,rhtxy,rho,rht,igq_fft)

      IF (slice)  STOP 'slice OK'
c-dw
c-t3e
      RETURN
      END SUBROUTINE cdngen
      END MODULE m_cdngen
