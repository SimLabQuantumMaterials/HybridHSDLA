      MODULE m_dimen7
      CONTAINS
      SUBROUTINE dimen7(
     X                  gw,gw_neigd,  
     <                  nop,k1d,k2d,k3d,n3d,n2d,
     <                  kq1d,kq2d,kq3d,kxc1d,kxc2d,kxc3d,
     <                  ntypd,natd,jmtd,ntypsd,nlhd,memd,
     <                  lmaxd,jspd,nvacd,
     <                  nvd,nv2d,nwdd,nkptd,nstd,
     <                  neigd,msh,ncvd,layerd,nlod,llod,
     <                  odd)

c
c This program reads the input files of the flapw-programm (inp & kpts)
c and creates a file 'fl7para' that contains dimensions 
c for the main flapw-programm.
c

      USE m_localsym
      USE m_socsym
      USE m_sssym
      USE m_spg2set
      USE m_constants, ONLY : pimach
      USE m_rwinp
      USE m_inpnoco
      USE m_julia
      USE m_od_kptsgen
      USE m_types,     ONLY : t_utype
      USE m_firstglance
      USE m_inv3
      USE m_rwsymfile
      USE m_strgndim
      USE m_convndim
      USE m_inpeigdim
      USE m_od_types, ONLY : od_dim

      IMPLICIT NONE

      INTEGER, INTENT (INOUT) :: gw,gw_neigd
c
c dimension-parameters for flapw:
c
      INTEGER, INTENT (OUT) :: jspd,nvacd
      INTEGER, INTENT (OUT) :: ntypd,natd
      INTEGER, INTENT (OUT) :: nop
      INTEGER, INTENT (OUT) :: n3d,n2d
      INTEGER, INTENT (OUT) :: nvd,nv2d
      INTEGER, INTENT (OUT) :: lmaxd
      INTEGER, INTENT (OUT) :: nkptd,nwdd
      INTEGER, INTENT (OUT) :: layerd
      INTEGER, INTENT (OUT) :: ntypsd,nlhd,memd
      INTEGER, INTENT (OUT) :: jmtd,ncvd,msh,nstd
      INTEGER, INTENT (OUT) :: neigd
      INTEGER, INTENT (OUT) :: k1d,k2d,k3d
      INTEGER, INTENT (OUT) :: kq1d,kq2d,kq3d
      INTEGER, INTENT (OUT) :: kxc1d,kxc2d,kxc3d
      INTEGER, INTENT (OUT) :: nlod,llod
      TYPE (od_dim), INTENT (OUT) :: odd
c
c-------------------------------------------------------------------
C ..  Local Scalars ..
      REAL thetad,xa,epsdisp,epsforce,rkmax,rmtmax,arltv1,arltv2,arltv3
      REAL dvac,dtild,scale,scpos,gmax,delgau,zc,tkb,alpha,spinf
      REAL e1s,e2s,chng,gmaxxc,theta,phi,omtil,tpi,tworkf,s,mix_b
      REAL e1_dos,e2_dos,sig_dos,r,d,thetaJ
      INTEGER ok,ntype,ilo,n,nstate,i,j,nstars,nstm,na,n1,jrc
      INTEGER isec1,ndir,lpr,nwd,lepr,kcrel,nkpt,nmop(3),nqpt,nmopq(3)
      INTEGER imix,layers,kk,nnne,maxiter,nop2,itmax,nopd,symfh
      INTEGER nmagn,nsh,mtypes
      INTEGER igrd,ndvgrd,idsprs
      CHARACTER(len=1) :: rw
      CHARACTER(len=3) :: latnam
      CHARACTER(len=4) :: namex,namgrp
      CHARACTER(len=7) :: symfn
      CHARACTER(len=12):: relcor
      LOGICAL l_u2f,l_f2u,l_bmt,l_soc,starcoeff,symor,l_kpts,l_qpts
      LOGICAL l_noco,l_ss,l_J,l_mperp,l_constr,l_inpexist,l_disp
      LOGICAL strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8
      LOGICAL form66,l_f,eonly,gauss,tria,frcor,slice,ctail,disp
      LOGICAL swsp,lflip,vacdos,integ,iplot,score,plpot,pallst,lwb
C ..
      REAL    :: amat(3,3),bmat(3,3),a1(3),a2(3),a3(3),locx(2),locy(2)
      REAL    :: qss(3),q(3)
      LOGICAL :: eig66(2),sso_opt(2) 

      CHARACTER(len=8) ::  name(10)
      CHARACTER(len=3), ALLOCATABLE :: noel(:)
      INTEGER, ALLOCATABLE :: lmax(:),mrot(:,:,:),jri(:)
      INTEGER, ALLOCATABLE :: ntypsy(:),neq(:),izlay(:,:)
      INTEGER, ALLOCATABLE :: ncst(:),lnonsph(:),nflip(:)
      INTEGER, ALLOCATABLE :: magtype(:),nmagtype(:)
      INTEGER, ALLOCATABLE :: nz(:),relax(:,:),llo(:,:),nlo(:),nlhtyp(:)
      INTEGER, ALLOCATABLE :: nlh(:),llh(:,:),nmem(:,:),mlh(:,:,:)
      LOGICAL, ALLOCATABLE :: l_geo(:),l_relax(:),error(:),l_magn(:)
      LOGICAL, ALLOCATABLE :: soc_opt(:)
      REAL,    ALLOCATABLE :: pos(:,:),rmt(:),zatom(:),dx(:),tau(:,:)
      REAL,    ALLOCATABLE :: ellow(:),elup(:),zelec(:),bmu(:),rkm(:)
      REAL,    ALLOCATABLE :: taual(:,:)
      REAL,    ALLOCATABLE :: alph(:),beta(:),b_con(:,:),M(:)
      COMPLEX, ALLOCATABLE :: clnu(:,:,:)

c+lda+u
      TYPE (t_utype), ALLOCATABLE :: lda_u(:)
c-lda+u
c-odim
      INTEGER ntp1,ii
      INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
c+odim

      EXTERNAL prp_xcfft_box,parawrite
C     ..
      tpi = 2 * pimach()

      DO i = 1,3
         a1(i) = 0
         a2(i) = 0
         a3(i) = 0
      ENDDO
c
c---> First, check whether an inp-file exists
c
      INQUIRE (file='inp',exist=l_inpexist)
      IF (.not.l_inpexist) THEN
!       INQUIRE (file='input',exist=l_inpexist)
!       IF (.not.l_inpexist) THEN
          STOP 'no inp- or input-file found!'
!       ENDIF
!       CALL inp_gen()
      ENDIF
c
c---> determine ntype,nop,natd,nwdd,nlod and layerd
c
      CALL first_glance(
     <                  ntype,nop,natd,nwdd,nlod,layerd,itmax,
     <                  l_kpts,l_qpts,nkpt,nmop,nqpt,nmopq)

      nlod = max(nlod,1)

      ALLOCATE (
     + lmax(ntype),ntypsy(natd),neq(ntype),nlhtyp(ntype),
     + rmt(ntype),zatom(ntype),jri(ntype),dx(ntype), 
     + nlo(ntype),llo(nlod,ntype),nflip(ntype),bmu(ntype),
     + noel(ntype),izlay(layerd,2),ncst(ntype),lnonsph(ntype),
     + taual(3,natd),pos(3,natd),rkm(nwdd),zelec(nwdd),
     + ellow(nwdd),elup(nwdd),nz(ntype),relax(3,ntype),
     + l_geo(ntype),soc_opt(ntype+2),alph(ntype),beta(ntype),
     + lda_u(ntype),l_relax(ntype),l_magn(ntype),M(ntype),
     + magtype(ntype),nmagtype(ntype),b_con(2,ntype),
     + clnu(1,1,1),nlh(1),llh(1,1),nmem(1,1),mlh(1,1,1),STAT=ok)
c
c---> read complete input and calculate nvacd,llod,lmaxd,jmtd,neigd and ncvd
c
      CALL rw_inp(
     >            6,ntype,natd,nwdd,layerd,nlod,5,'R',
     < dvac,dtild,scale,scpos,gmax,rkm,delgau,zc,tkb,alpha,spinf,
     < e1s,e2s,isec1,ndir,jspd,lpr,nwd,lepr,kcrel,itmax,imix,layers,
     < l_u2f,l_f2u,l_bmt,kk,nnne,maxiter,latnam,noel,namex,namgrp,
     < relcor,
     < strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8,form66,
     < l_f,eonly,eig66,gauss,tria,frcor,slice,ctail,disp,swsp,lflip,
     < vacdos,integ,iplot,score,plpot,pallst,a1,a2,a3,rmt,dx,ellow,
     < elup,zelec,bmu,ncst,jri,lmax,lnonsph,nflip,izlay,
     < name,igrd,ndvgrd,idsprs,lwb,chng,gmaxxc,l_soc,soc_opt,theta,phi,
     < taual,
     < ntype,neq,nz,xa,relax,thetad,epsdisp,epsforce,nlo,llo,tworkf,
     < nstars,nstm,starcoeff,locx,locy,l_noco,l_J,l_geo,e1_dos,e2_dos,
     < sig_dos,lda_u,gw,gw_neigd,odd)

c---> pk non-collinear
c---> read the angle and spin-spiral information from nocoinp
      qss(1) = 0.0
      qss(2) = 0.0
      qss(3) = 0.0
      l_ss = .false.
      IF (l_noco) THEN 
        IF (.not.l_J) THEN
         CALL inpnoco(
     >    natd,24,itmax,isec1,secvar,gauss,tria,
     >    ctail,score,ntype,neq,taual,nstm,starcoeff,l_f,l_J,l_soc,
     <    l_ss,l_mperp,l_constr,mix_b,alph,beta,b_con,
     <    l_relax,qss,thetaJ,l_magn,nmagn,M,mtypes,
     <    magtype,nmagtype,nsh,l_disp,sso_opt)
        ELSE
         l_ss= .true.
        ENDIF
      ENDIF

      nvacd = 2
      IF (zrfs .OR. invs .OR. odd%d1) nvacd = 1
      llod = 0
      lmaxd = 0
      jmtd = 0
      rmtmax = 0.0
      neigd = 0
      nstd = 0
      DO n = 1,ntype
        nstd = max(ncst(n),nstd)
        lmaxd = max(lmax(n),lmaxd)
        jmtd = max(jri(n),jmtd)
        rmtmax =  max(rmtmax,rmt(n))
        DO ilo = 1,nlo(n)
c+apw
          IF (llo(ilo,n).LT.0) THEN
             llo(ilo,n) = -llo(ilo,n) - 1
#ifndef CPP_APW
             llo(ilo,n) = mod(llo(ilo,n),10)
#endif
          ELSE
             neigd = neigd + neq(n)*(2*abs(llo(ilo,n)) +1)
          ENDIF
c-apw
          llod = max(abs(llo(ilo,n)),llod)
        ENDDO
        nstate = 4
        IF ((nz(n).GE.21.AND.nz(n).LE.29) .OR. 
     +      (nz(n).GE.39.AND.nz(n).LE.47) .OR.
     +      (nz(n).GE.57.AND.nz(n).LE.79)) nstate = 9
        IF ((nz(n).GE.58.AND.nz(n).LE.71) .OR.
     +      (nz(n).GE.90.AND.nz(n).LE.103)) nstate = 16
        neigd = neigd + nstate*neq(n)
c
      ENDDO
      
!      IF (mod(lmaxd,2).NE.0) lmaxd = lmaxd + 1
      IF (2*neigd.LT.zc) THEN
        WRITE(6,*) neigd,' states estimated in dimen7 ...'
        neigd = NINT(0.75*zc)
        WRITE(6,*) 'changed neigd to ',neigd
      ENDIF
      IF (l_soc .and. (.not. l_noco)) neigd=2*neigd 
      IF (l_soc .and. l_ss) neigd=(3*neigd)/2  
       ! not as accurate, but saves much time

      rmtmax = rmtmax*gmax
      CALL convn_dim(
     >               rmtmax,
     <               ncvd)
!
! determine core mesh
!
      msh = 0
      DO n = 1,ntype
         r = rmt(n)
         d = exp(dx(n))
         jrc = jri(n)
         DO WHILE (r < rmt(n) + 20.0)
            jrc = jrc + 1
            r = r*d
         ENDDO
         msh = max( msh, jrc ) 
      ENDDO
c
c ---> now, set the lattice harmonics, determine nlhd
c
      DO i = 1,3
         amat(i,1) = a1(i)*scale
         amat(i,2) = a2(i)*scale
         amat(i,3) = a3(i)*scale
      ENDDO
      CALL inv3(amat,bmat,omtil)
      IF (film) omtil = omtil/amat(3,3)*dvac
c-odim
      IF (odd%d1) omtil = amat(3,3)*pimach()*(dvac**2)/4.
c+odim
      DO j = 1,3
         DO i = 1,3
            bmat(i,j) = tpi*bmat(i,j)
         ENDDO
      ENDDO

      na = 0
      DO n = 1,ntype
        DO n1 = 1,neq(n)
            na = na + 1
            IF (film) taual(3,na) = taual(3,na)/a3(3)
            DO i = 1,3
              s = 0
              DO j = 1,3
                 s = s + amat(i,j)*taual(j,na)
              ENDDO
              pos(i,na) = s
           ENDDO 
        ENDDO
        zatom(n) = real( nz(n) )
      ENDDO
      ALLOCATE (mrot(3,3,nop),tau(3,nop))
      IF (namgrp.EQ.'any ') THEN
         nopd = nop ; rw = 'R'
         symfh = 94 ; symfn = 'sym.out'
         CALL rw_symfile(
     >                   rw,symfh,symfn,nopd,bmat,
     X                   mrot,tau,nop,nop2,symor)
      ELSE
         CALL spg2set(
     >                nop,zrfs,invs,namgrp,latnam,
     <                mrot,tau,nop2,symor)
      ENDIF
      ntypsd = 0
      IF (.NOT.odd%d1) THEN
        CALL local_sym(
     >               lmaxd,lmax,nop,mrot,tau,
     >               natd,ntype,neq,amat,bmat,taual,
     <               nlhd,memd,ntypsd,.true.,
     =               nlhtyp,ntypsy,nlh,llh,nmem,mlh,clnu)
c-odim
      ELSEIF (odd%d1) THEN
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
     >      lmaxd,lmx1,nop,mrot,tau,
     >      natd,ntp1,nq1,amat,bmat,taual,
     <      nlhd,memd,ntypsd,.true.,
     =      nlhtp1,ntypsy,nlh,llh,nmem,mlh,clnu)
        ii = 1
        DO i = 1,ntype
          nlhtyp(i) = nlhtp1(ii)
          ii = ii + neq(i)
        END DO
        DEALLOCATE (nq1,lmx1,nlhtp1)
      END IF
c+odim
!
! Check if symmetry is compatible with SOC or SSDW
!
      IF (l_soc .and. (.not.l_noco)) THEN  
        ! test symmetry for spin-orbit coupling
        ALLOCATE ( error(nop) )
        CALL soc_sym(
     >               nop,mrot,theta,phi,amat,
     <               error)
        IF ( ANY(error(:)) ) THEN
          WRITE(*,fmt='(x)')                                    
          WRITE(*,fmt='(A)')                                    
     &     'Symmetry incompatible with SOC spin-quantization axis ,'  
          WRITE(*,fmt='(A)')                                   
     &     'do not perform self-consistent calculations !'    
          WRITE(*,fmt='(x)')                                 
          IF ( eonly .or. (l_soc.and.l_ss) ) THEN  ! .or. ... .or. ...  
            CONTINUE 
          ELSE 
            IF (itmax>1) THEN
              STOP 'dimen7 : symmetry & SOC'
            ENDIF 
          ENDIF 
        ENDIF           
        DEALLOCATE ( error )
      ENDIF
c--- J<
      IF(.not.l_J) THEN
c--- J>
      IF (l_ss) THEN  ! test symmetry for spin-spiral
        ALLOCATE ( error(nop) )
        CALL ss_sym(
     >              nop,mrot,qss,
     <              error)
        IF ( ANY(error(:)) ) STOP 'dimen7 : symmetry & SSDW'
        DEALLOCATE ( error )
      ENDIF
c--- J<
      ENDIF
c--- J>
!
! Dimensioning of the stars
!
      IF (film.OR.(namgrp.ne.'any ')) THEN
         CALL strgn1_dim(
     >                gmax,bmat,invs,zrfs,mrot,tau,nop,nop2,tpi,
     <                k1d,k2d,k3d,n3d,n2d,odd)

      ELSE
         CALL strgn2_dim(
     >                gmax,bmat,invs,zrfs,mrot,tau,nop,tpi,
     <                k1d,k2d,k3d,n3d,n2d)
         odd%n2d = n2d
         odd%nq2 = n2d
         odd%nop = nop
      ENDIF

      IF ( gmaxxc .le. 10.0**(-6) ) THEN
         WRITE (6,'(" gmaxxc=0 : gmaxxc=gmax choosen as default",
     +              " value")')
         WRITE (6,'(" concerning memory, you may want to choose",
     +              " a smaller value for gmax")')
         gmaxxc=gmax
      END IF

      CALL prp_xcfft_box(gmaxxc,bmat,kxc1d,kxc2d,kxc3d)
c
c k-point generator provides kpts-file, if it's missing:
c
      IF (.not.l_kpts) THEN
       IF (.NOT.odd%d1) THEN
         IF (l_J) THEN
         CALL julia(
     >              1,   1,mrot,tau,amat,bmat,omtil,jspd,latnam,
     >              film,tria,invs,invs2,l_ss,l_soc,ndir,name,
     X              nkpt,nmop,.false.)
         ELSE
         CALL julia(
     >              nop,nop2,mrot,tau,amat,bmat,omtil,jspd,latnam,
     >              film,tria,invs,invs2,l_ss,l_soc,ndir,name,
     X              nkpt,nmop,.false.)
         ENDIF
       ELSE
        CALL od_kptsgen (nkpt)
       ENDIF
      ELSE
        IF(gw.eq.2) THEN
          neigd=gw_neigd ! gw_neigd is ignored if gw=0 or gw=1.
          INQUIRE(file='QGpsi',exist=l_kpts) ! Use QGpsi if it exists otherwise use kpts-file.
          IF(l_kpts) THEN
            WRITE(6,*)
     &        'QGpsi exists and will be used to generate kpts-file'
            OPEN (15,file='QGpsi',form='unformatted',status='old',
     &        action='read')
            OPEN (41,file='kpts',form='formatted',status='unknown')
            REWIND(41)
            READ (15) nkpt
            WRITE (41,'(i5,f20.10)') nkpt,1.0
            DO n = 1, nkpt
              READ (15) q
              WRITE (41,'(4f10.5)') MATMUL(TRANSPOSE(amat),q)/scale,1.0
              READ (15)
            ENDDO
            CLOSE (15)
            CLOSE (41)
          ENDIF
        ENDIF
      ENDIF
c
c Using the k-point generator also for creation of q-points for the
c J-constants calculation:
      IF(.not.l_qpts)THEN
        CALL julia(
     >             nop,nop2,mrot,tau,amat,bmat,omtil,jspd,latnam,
     >             film,tria,invs,invs2,.false.,.false.,ndir,name,
     X             nqpt,nmopq,.true.)
      ENDIF

c
c now proceed as usual
c
      rkmax = max( rkm(1),rkm(nwdd) )
      CALL inpeig_dim(
     >                rkmax,nwdd,latnam,amat,bmat,film,l_ss,qss,odd,l_J,
     <                nkptd,nvd,nv2d,kq1d,kq2d,kq3d)
      layerd = max(layerd,1)
      nstd = max(nstd,30)
      ntypd = ntype
      IF (l_noco) neigd = 2*neigd

      nlod = max(nlod,2) ! for chkmt
      CALL parawrite(
     >               nop,k1d,k2d,k3d,n3d,n2d,nlod,llod,
     >               kq1d,kq2d,kq3d,kxc1d,kxc2d,kxc3d,
     >               ntypd,natd,jmtd,ntypsd,nlhd,memd,lmaxd,
     >               jspd,nvacd,nvd,nv2d,nwdd,msh,
     >               nstd,nkptd,neigd,ncvd,layerd,odd)

c
      DEALLOCATE( 
     + lmax,ntypsy,neq,nlhtyp,rmt,zatom,jri,dx,nlo,llo,nflip,bmu,noel,
     + izlay,ncst,lnonsph,taual,pos,rkm,zelec,ellow,elup,nz,relax,
     + l_geo,soc_opt,alph,beta,lda_u,l_relax,l_magn,M,b_con,clnu,nlh,
     + llh,nmem,mlh,magtype,nmagtype)
c
      RETURN
      END SUBROUTINE dimen7
      END MODULE m_dimen7
