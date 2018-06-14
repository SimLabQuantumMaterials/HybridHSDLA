      MODULE m_vgen
      CONTAINS
      SUBROUTINE vgen(
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
     >                invtab,phi2,odi,ods,odl,odg,plpot,irank,isize,
     <                vr0,te_vcoul,te_veff,te_exc,l_noco,l_ss,qss)
c     ***********************************************************
c     FLAPW7 film potential generator                           *
c     composed July 1986, revised and tested May 1987           *
c     spin-polarized added may 1987                             *
c                         c.l.fu,  R. Podloucky                 *
c     ***********************************************************
c     calculates the density-potential integrals needed for the
c     total energy
c     TE_VCOUL  :   charge density-coulomb potential integral
c     TE_VEFF:   charge density-effective potential integral
c     TE_EXC :   charge density-ex-corr.energy density integral
c     ***********************************************************
c     VR0     average coulomb potential contributes to the Madelung term
c     R.Pentcheva 09.05.96
c     ***********************************************************
c
#include"cpp_double.h"
      USE m_constants, ONLY : pimach
      USE m_vmts
      USE m_intnv
      USE m_vmtxcg
      USE m_vvacxc
      USE m_vvacxcg
      USE m_visxc
      USE m_visxcg
      USE m_vvac
      USE m_vvacis
      USE m_vvacxy
      USE m_vintcz
      USE m_checkdop
      USE m_loddop
      USE m_wrtdop
      USE m_qfix
      USE m_od_types
      USE m_od_vvac
      USE m_od_vvacis
      USE m_cylpts
      USE m_cputime
      USE m_outtime
      USE m_convol
      USE m_xyavden
      USE m_psqpw
      USE m_potmod
!pg      USE m_intgr, ONLY : intgr3

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd
      INTEGER, INTENT (IN) :: ncvd,nmzd,k1d,k2d,k3d,kxc1d,kxc2d,kxc3d
      INTEGER, INTENT (IN) :: nop,n2d,natd,nn2d,nmzxyd,nspd
      INTEGER, INTENT (IN) :: nq2,nq3,nmzxy,nmz,mx3,nk2,nk1,jspins,ntype
      INTEGER, INTENT (IN) :: nvac,icorr,nsymt,nop2,it,nwd
      INTEGER, INTENT (IN) :: kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft
      INTEGER, INTENT (IN) :: kmxxc_fft,kimax,kimax2,irank,isize
      INTEGER, INTENT (IN) :: idsprs0,idsprsl,idsprsi,idsprsv
      INTEGER, INTENT (IN) :: igrd,ndvgrd,isprsv,krla
      REAL,    INTENT (IN) :: chng,sprsv,sigma,zsigma
      REAL,    INTENT (IN) :: z1,area,dvac,delz,omtil,volint,vol
      REAL,    INTENT (OUT):: te_vcoul,te_veff,te_exc
      LOGICAL, INTENT (IN) :: reap,lwb,invs2,invs,zrfs,total,film,vchk
      LOGICAL, INTENT (IN) :: symor,l_noco,l_ss,plpot
      INTEGER, INTENT (INOUT) :: idsprs
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),lmax(ntypd),jri(ntypd)
      INTEGER, INTENT (IN) :: nlh(ntypsd),neq(ntypd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),kv2(2,n2d),kv3(3,n3d)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),ncv(ntypd),ngopr(natd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      INTEGER, INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),qss(3)
      REAL,    INTENT (IN) :: pos(3,natd),taual(3,natd)
      REAL,    INTENT (IN) :: volmts(ntypd),rmsh(jmtd,ntypd),rmt(ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),sk2(n2d),sk3(n3d),zatom(ntypd)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL,    INTENT (IN) :: tau(3,nop),sig_b(2)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1),pgft2xy(0:nn2d-1)
      REAL,    INTENT (IN) :: pgft2x(0:nn2d-1),pgft2y(0:nn2d-1)
      REAL,    INTENT (IN) :: pgft2xx(0:nn2d-1),pgft2yy(0:nn2d-1)
      REAL,        INTENT (OUT) :: vr0(ntypd)
      CHARACTER*2, INTENT (IN)  :: namat(0:103)
c-odim
      REAL,    INTENT (IN) :: phi2(n2d)
      TYPE (od_inp), INTENT (INOUT) :: odi
      TYPE (od_sym), INTENT (INOUT) :: ods
      TYPE (od_lda), INTENT (INOUT) :: odl
      TYPE (od_gga), INTENT (INOUT) :: odg
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX vintcza,xint,rhobar
      INTEGER i,i3,irec2,irec3,ivac,j,js,k,k3,lh,n,nzst1
      INTEGER imz,imzxy,i2d,ichsmrg,ivfft,nt,npd,nat
      INTEGER ifftd,ifftd2, ifftxc3d,iter,datend
      INTEGER i_sm,n_sm,i_sta,i_end
      REAL ani,g3,signum,time1,time2,z,rhmn,fix,fpi,sfp,mfie
      REAL sig1dh,vz1dh,zat_l(ntypd),delta,deltb,corr
      CHARACTER*8 dop,iop,name(10)
      LOGICAL l_pottot
      LOGICAL, PARAMETER :: l_xyav=.false.
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: alphm(:,:)
      COMPLEX, ALLOCATABLE :: excpw(:),excxy(:,:,:),vpw_w(:,:),psq(:)
      REAL,    ALLOCATABLE :: vbar(:),af1(:),bf1(:),xp(:,:)
c.....density
      REAL,    ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
c ff
      COMPLEX, ALLOCATABLE :: cdom(:)
      COMPLEX, ALLOCATABLE :: cdomvz(:,:),cdomvxy(:,:,:)
c.....potential
      REAL,    ALLOCATABLE :: vr(:,:,:,:),vz(:,:,:)
      COMPLEX, ALLOCATABLE :: vpw(:,:),vxy(:,:,:,:)
c.....energy density
      REAL,    ALLOCATABLE :: excz(:,:),excr(:,:,:)
!
! if you want to calculate potential gradients
!
!pg      COMPLEX vlm(-lmaxd:lmaxd,0:lmaxd,ntypd)
!pg      REAL fint,f(jmtd)
!pg      INTEGER ns, nl, l, jm, m, mb
C     ..
C     .. External Subroutines ..
      EXTERNAL cfft,points,sphpts,vmtxc
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aimag,cmplx,real,sqrt
C     ..
c     ----> note the following conventions:
c     ivac=2: lower (negative z) vacuum
c     ivac=1: upper (positive z) vacuum
c     units: hartrees
c
      ALLOCATE ( alphm(n2d,2),excpw(n3d),excxy(nmzxyd,odi%n2d-1,2),
     +           vbar(jspd),af1(3*k3d),bf1(3*k3d),xp(3,nspd),
     +           rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd),
     +           qpw(n3d,jspd),rhtxy(nmzxyd,odi%n2d-1,2,jspd),
     +           vr(jmtd,0:nlhd,ntypd,jspd),vz(nmzd,2,jspd),
     +           vpw(n3d,jspd),vxy(nmzxyd,odi%n2d-1,2,jspd),
     +           excz(nmzd,2),excr(jmtd,0:nlhd,ntypd),
     +           vpw_w(n3d,jspd),psq(n3d) )
      IF (l_noco) THEN
        ALLOCATE ( cdom(n3d),
     +             cdomvz(nmzd,2),cdomvxy(nmzxyd,odi%n2d-1,2) )
      ELSE
        ALLOCATE ( cdom(1),cdomvz(1,1),cdomvxy(1,1,1) )
      END IF
c
      fpi = 4 * pimach()
      sfp = sqrt( fpi )
      IF (irank == 0) THEN
c ff
      IF (l_noco) THEN
        ! nt = 70
        ! OPEN (nt,file='cdn',form='unformatted',status='old')
        ! In the non-collinear case |m| is calculated from mx,my,mz
        ! in "visxc(g)","vvacxc(g)" instead of using the "cdn"-file.
        ! It is done this way as accuracy is lost when transforming
        ! |m| back to g-space in "rhodirgen" .
        nt = 26
        OPEN (nt,file='rhomat_inp',form='unformatted',status='old')
      ELSE
         nt = 71
         OPEN (nt,file='cdn1',form='unformatted',status='old')
      ENDIF
c
c --  total = .false. and reap = .false. means, that we now calculate
c     the potential from the output density
c
      IF ((.not.total).AND.(.not.reap)) THEN
         IF (l_noco) THEN
           STOP 'vgen:1'
         ENDIF
         CALL loddop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >               jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >               nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <               iop,dop,iter,rho,qpw,rht,rhtxy,name)
      ENDIF
c
c --> load density and fix
c
      CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
      IF (l_noco) THEN
        READ (nt,iostat=datend) (cdom(k),k=1,nq3)
        IF (datend == 0) THEN
          IF (film) THEN
            READ (nt) ((cdomvz(i,ivac),i=1,nmz),ivac=1,nvac)
            READ (nt) (((cdomvxy(i,j-1,ivac),i=1,nmzxy),j=2,odi%nq2),
     +                                               ivac=1,nvac)
          END IF
        ELSE
        ! (datend < 0)  =>  no off-diagonal magnetisation stored
        !                   in "rhomat_inp"
          IF (datend > 0) THEN
            STOP 'vgen: error reading "rhomat_inp" '
          END IF
          DO k=1,nq3
            cdom(k)= CMPLX(0.,0.)
          END DO
          IF (film) THEN
            DO ivac=1,nvac
              DO i=1,nmz
                cdomvz(i,ivac)= CMPLX(0.,0.)
              END DO
              DO i=1,nmzxy
                DO j=2,odi%nq2
                  cdomvxy(i,j-1,ivac)= CMPLX(0.,0.)
                END DO
              END DO
            END DO
          END IF
        END IF
      END IF

      IF (.not.l_xyav) THEN
      CALL qfix(
     >          k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >          nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,nmzxy,n2d,
     >          ntype,neq,volmts,taual,z1,vol,volint,nq2,invtab,
     >          symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >          nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     >          qpw,rhtxy,rho,rht,odi,
     <          fix)
      ENDIF
      IF (total.OR.reap) REWIND nt
      CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt,
     >            iop,dop,iter,rho,qpw,rht,rhtxy,name)
      IF (l_noco) THEN
        WRITE (nt) (cdom(k),k=1,nq3)
        IF (film) THEN
          WRITE (nt) ((cdomvz(i,ivac),i=1,nmz),ivac=1,nvac)
          WRITE (nt) (((cdomvxy(i,j-1,ivac),i=1,nmzxy),j=2,odi%nq2),
     +                                              ivac=1,nvac)
        ENDIF
      ENDIF
c
      WRITE (6,FMT=8000)
 8000 FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)
      DO 20 js = 1,jspins
         DO 10 k = 1,n3d
            vpw(k,js) = cmplx(0.,0.)
   10    CONTINUE
   20 CONTINUE
c     ---> perform spin summation of charge densities
c     ---> for the calculation of the coulomb potentials
      IF ( (jspins.EQ.2).AND.(.not.l_xyav) ) THEN
         nat = 1
         DO 50 n = 1,ntype
            DO 40 lh = 0,nlh(ntypsy(nat))
               DO 30 i = 1,jri(n)
                  rho(i,lh,n,1) = rho(i,lh,n,1) + rho(i,lh,n,jspins)
   30          CONTINUE
   40       CONTINUE
            nat = nat + neq(n)
   50    CONTINUE
         DO 60 i = 1,nq3
            qpw(i,1) = qpw(i,1) + qpw(i,jspins)
   60    CONTINUE
         IF (film) THEN
            DO 100 ivac = 1,nvac
               DO 80 j = 1,odi%nq2 - 1
                  DO 70 i = 1,nmzxy
                     rhtxy(i,j,ivac,1) = rhtxy(i,j,ivac,1) +
     +                                   rhtxy(i,j,ivac,jspins)
   70             CONTINUE
   80          CONTINUE
               DO 90 i = 1,nmz
                  rht(i,ivac,1) = rht(i,ivac,1) + rht(i,ivac,jspins)
   90          CONTINUE
  100       CONTINUE
         END IF
      END IF
c
c     ************** coulomb potential ***********************
      IF (l_xyav) THEN
!        rho(:,:,:,1) = rho(:,:,:,jspins)     ! just if you want to take 
!        qpw(:,1) = qpw(:,jspins)             ! the other spin direction
!        rhtxy(:,:,:,1) = rhtxy(:,:,:,jspins)
!        rht(:,:,1) = rht(:,:,jspins)
        zat_l(:) = 0.          ! for xy-averaged densities do not
      ELSE                     ! include nuclear charge
        zat_l(:) = zatom(:)
      ENDIF

c     ----> create pesudo-charge density coefficients
      CALL cpu_time(time1)
      ENDIF ! (irank == 0)
      CALL psqpw(irank,isize,
     >           ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,memd,nmzd,
     >           ntype,nq3,jri,ntypsy,lmax,jspd,ncvd,
     >           rmt,rmsh,dx,zat_l,volmts,sfp,sk3,nstr,
     >           omtil,dvac,z1,area,film,nvac,neq,nmz,delz,
     >           clnu,mlh,nmem,llh,nlh,bmat,kv3,ig2,
     >           nop,natd,symor,mrot,tau,taual,invtab,
     >           ncv,qpw,rho,rht,sigma,vol,odi,ods,l_xyav,
     <           psq)
      IF (irank == 0) THEN
      CALL cpu_time(time2)
      CALL outtime('psqpw:', time2-time1)
#ifdef CPP_GF
      OPEN (12,file='gf_pseudo',form='unformatted')
      WRITE (12) n3d
      WRITE (12) psq
      CLOSE (12)
#endif

      IF (l_xyav) THEN        ! write out xy-averaged density & stop
        CALL xy_av_den(
     >                 n3d,k3d,nq3,nmzd,nmz,dvac,delz,
     >                 area,ig2,kv3,amat,psq,rht(1,1,1))
        STOP 'xy-averaged density calculated'
      ENDIF

c     ------------------------------------------
      IF (odi%d1) THEN
c-odim
         CALL cpu_time (time1)

c---> generates the m=0,gz=0 component of the vacuum potential
         CALL od_vvac(
     >        n3d,k1d,k2d,k3d,n2d,jspd,nvac,
     >        ig,ig2,sk2,phi2,nq3,nq2,
     >        z1,nmz,nmzd,delz,
     >        psq,rht,
     <        vz)

c---> generation of the vacuum warped potential components and
c---> interstitial pw potential
c---> vvacxy_5.F is a symmetrized potential generator

         CALL od_vvacis(
     >        odi%n2d,jspd,nmzd,nmzxyd,nmzxy,odi%nq2,
     >        odi%kv,nvac,z1,delz,odi%M,nq3,nstr2,
     >        k3d,bmat,sk3,amat,nq2,odi%nst2,
     >        k1d,k2d,sk2,phi2,odi%ig,ig,ig2,n2d,n3d,
     >        igfft2,pgfft2,kimax2,kv2,kv3,omtil,
     >        rht,rhtxy,psq,vz,zrfs,invs,
     >        kimax,nstr,igfft,pgfft,
     <        vxy,vpw)

         CALL cpu_time (time2)

         CALL outtime('coulomb potential:',
     +        time2-time1)

c+odim
      ELSEIF (film .AND. .NOT.odi%d1) THEN
c     ----> potential in the  vacuum  region
c        
         CALL cpu_time(time1)
         CALL vvac(
     >             nmzd,jspd,n3d,nmz,nvac,dvac,
     >             z1,nq3,zrfs,invs,delz,nstr,bmat,ig2,kv3,
     >             psq,rht,sigma,zsigma,sig_b,area,
     <             vz,rhobar,sig1dh,vz1dh)
         CALL vvacis(
     >               k1d,k2d,k3d,n3d,n2d,jspd,nmzxyd,nvac,nmzxy,nq2,
     >               zrfs,mx3,z1,delz,bmat,nstr2,nstr,ig,sk2,kv2,
     >               psq,
     <               vxy)
         CALL vvacxy(
     >               n2d,jspd,nmzxyd,nmzxy,nq2,nvac,z1,delz,invs,sk2,
     >               rhtxy,
     X               vxy,
     <               alphm)
         CALL cpu_time(time2)
         CALL outtime('p vac', time2-time1)
      END IF
c     ------------------------------------------
c     ----> potential in the  interstitial  region
      WRITE (6,FMT=8010)
 8010 FORMAT (/,5x,'coulomb potential in the interstitial region:')
      IF (film .AND. .NOT.odi%d1) THEN
c           -----> create v(z) for each 2-d reciprocal vector
         ivfft =  3*k3d 
!         ivfft = 2*mx3 - 1
         ani = 1.0/real(ivfft)
         DO 130 irec2 = 1,nq2
            i = 0
            DO 110 i3 = 0,ivfft - 1
               i = i + 1
               z = amat(3,3)*i3*ani
               IF (z.GT.amat(3,3)/2.) z = z - amat(3,3)
               vintcza = vintcz(
     >                          n3d,nmzxyd,n2d,jspd,k1d,k2d,k3d,nmzd,
     >                          z,irec2,sk3,sk2,nstr2,kv2,nstr,bmat,
     >                          z1,invs,zrfs,delz,nmzxy,mx3,ig,
     >                          psq,vxy,vz,rhobar,sig1dh,vz1dh,alphm)
               af1(i) = real(vintcza)
               bf1(i) = aimag(vintcza)
  110       CONTINUE
!            IF (irec2 == 1) THEN                 ! smooth potential
!              DO i_sm = 1, ivfft
!                z = (i_sm-1)*ani
!                IF (z > 0.5) THEN
!                  i_sta = i_sm -1
!                  i_end = i_sm 
!                  EXIT
!                ENDIF
!              ENDDO
!              delta = af1(i_end) - af1(i_sta) 
!              deltb = bf1(i_end) - bf1(i_sta) 
!              delta = delta  * ivfft / (ivfft - 1)
!              deltb = deltb  * ivfft / (ivfft - 1)
!              DO i_sm = 1, ivfft
!                z = (i_sm-1)*ani
!                IF (z > 0.5) z = z - 1.0
!                af1(i_sm) = af1(i_sm) + z * delta
!                bf1(i_sm) = bf1(i_sm) + z * deltb
!              ENDDO
!            ENDIF
c        --> 1-d fourier transform and store the coefficients in vpw( ,1)
            CALL cfft(af1,bf1,ivfft,ivfft,ivfft,-1)
!            delta = ivfft * delta * 2 / fpi ! * amat(3,3)**2 * ani
            i = 0
            DO 120 i3 = 0,ivfft - 1
               k3 = i3
               IF (k3 > FLOOR(ivfft/2.) ) k3 = k3 - ivfft
               i = i + 1
               IF ((k3.GE.-mx3).AND.(k3.LE.mx3)) THEN
                 irec3 = ig(kv2(1,irec2),kv2(2,irec2),k3)

!                 IF ( (irec2 == 1).AND.(i3 > 0) ) THEN                 ! smooth potential
!                   corr = 2.0*mod(abs(k3),2) - 1.0
!                   bf1(i) = bf1(i) + delta * corr / k3
!                 ENDIF

c       ----> only stars within g_max sphere (shz oct.97)
                 IF (irec3.ne.0) THEN
c
                    xint = cmplx(af1(i),bf1(i))*ani
                    nzst1 = nstr(irec3)/nstr2(irec2)
                    vpw(irec3,1) = vpw(irec3,1) + xint/nzst1
                 END IF
               ENDIF
  120       CONTINUE
  130    CONTINUE
      ELSEIF (.NOT.film) THEN
         vpw(1,1) = cmplx(0.0,0.0)
         DO 140 irec3 = 2,nq3
            g3 = sk3(irec3)
            vpw(irec3,1) = fpi*psq(irec3)/ (g3*g3)
  140    CONTINUE
      END IF

#ifdef CPP_GF
      OPEN  (12,file='gf_coulpot',form='unformatted')
      WRITE (12) n3d
      WRITE (12) vpw
      CLOSE (12)
#endif
      CALL cpu_time(time1)
      CALL outtime('p int', time1-time2)
      ENDIF ! irank == 0
c     --------------------------------------------
c     ---> potential in the muffin-tin spheres

      CALL vmts(irank,isize,
     >          n3d,jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,
     >          ntype,nq3,zatom,dx,rmsh,jri,rmt,sk3,lmax,
     >          clnu,mlh,nmem,llh,nstr,ntypsy,nlh,invtab,
     >          nop,natd,neq,kv3,mrot,bmat,tau,taual,symor,
     >          vpw,rho,odi,ods,
     <          vr)
c     --------------------------------------------
      IF (irank == 0) THEN
      CALL cpu_time(time2)
      CALL outtime('p mts', time2-time1)
c     ---> check continuity of coulomb potential

      IF (vchk) THEN
c           ----> vacuum boundaries
         IF (film .AND. .NOT.odi%d1) THEN
            npd = min(nspd,25)
            CALL points(xp,npd)
            DO 160 ivac = 1,nvac
              signum = 3. - 2.*ivac
              DO 150 j = 1,npd
                 xp(3,j) = signum*z1/amat(3,3)
  150         CONTINUE
              CALL checkdop(
     >                    xp,npd,0,0,ivac,1,1,.false.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    vpw,vr,vxy,vz,odi,ods)
  160       CONTINUE
         ELSEIF (odi%d1) THEN
c-odim
            npd = min(nspd,25)
            CALL cylpts(xp,npd,z1)
c           DO j = 1,npd
c              xp(1,j) = xp(1,j)/amat(1,1)
c              xp(2,j) = xp(2,j)/amat(2,2)
c           ENDDO
            CALL checkdop(
     >           xp,npd,0,0,nvac,1,1,.false.,nspd,jmtd,
     >           memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >           jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >           nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >           tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >           vpw,vr,vxy,vz,odi,ods)
c+odim
         END IF
c           ----> m.t. boundaries
         nat = 1
         DO 170 n = 1,ntype
            CALL sphpts(xp,nspd,rmt(n),pos(1,nat))
            CALL checkdop(
     >                    xp,nspd,n,nat,0,-1,1,.false.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    vpw,vr,vxy,vz,odi,ods)
            nat = nat + neq(n)
  170    CONTINUE
      END IF
C
C========TOTAL==============================================
C
!      IF (l_xyav) THEN        ! write out xy-averaged potential & stop
!        CALL xy_av_den(
!     >                 n3d,k3d,nq3,nmzd,nmz,dvac,delz,
!     >                 area,ig2,kv3,amat,vpw,vz(1,1,1))
!        STOP 'xy-averaged potential calculated'
!      ENDIF

      IF (total) THEN

c
c      ---> AVERAGE COULOMB POTENTIAL ON THE SPHERE 
C          FOR CALCULATING THE MADELUNG TERM in totale.f
c           r=Rmt
         DO 175 n = 1,ntype
            vr0(n) = vr(jri(n),0,n,1)
 175     CONTINUE
C
C     CALCULATE THE INTEGRAL OF n*Vcoulomb
C
         WRITE (6,FMT=8020)
         WRITE (16,FMT=8020)
 8020    FORMAT (/,10x,'density-coulomb potential integrals',/)
c
c       interstitial first
c
c       convolute ufft and pot: F(G) = \sum_(G') U(G - G') V(G')
c
        CALL convol(
     >              k1d,k2d,k3d,n3d,
     <              vpw_w,
     >              vpw,nq3,
     =              kimax,igfft,pgfft,ufft)
c
        IF (jspins.EQ.2) CALL CPP_BLAS_ccopy(nq3,vpw_w(1,1),1,
     +                                       vpw_w(1,jspins),1)
c
        te_vcoul = 0.0
        CALL int_nv(
     >             n3d,n2d,nmzxyd,nmzd,jmtd,ntypd,nlhd,ntypsd,
     >             qpw,vpw_w,
     >             rhtxy,vxy,
     >             rht,vz,
     >             rho,vr,
     >             area,delz,omtil,rmsh,dx,
     >             jri,neq,nlh,nstr2,ntypsy,natd,
     >             nq3,
     >             nmz,nmzxy,nq2,ntype,nvac,
     >             zrfs,film,invs,z1,odi,
     <             te_vcoul)

      WRITE (6,FMT=8030) te_vcoul
      WRITE (16,FMT=8030) te_vcoul
 8030 FORMAT (/,10x,'total density-coulomb potential integral :',
     +     t40,f20.10)

      END IF
C
C==========END TOTAL===================================================
c
c     ----> reload the density for calculating vxc (for spin-pol. case)
c
      IF (jspins.EQ.2) THEN
c
         REWIND nt
         CALL loddop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >               jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >               nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <               iop,dop,iter,rho,qpw,rht,rhtxy,name)
         CLOSE (nt)

         nat = 1
         DO 200 n = 1,ntype
            DO 190 lh = 0,nlh(ntypsy(nat))
               DO 180 i = 1,jri(n)
                  vr(i,lh,n,jspins) = vr(i,lh,n,1)
  180          CONTINUE
  190       CONTINUE
            nat = nat + neq(n)
  200    CONTINUE
         DO 210 i = 1,nq3
            vpw(i,jspins) = vpw(i,1)
  210    CONTINUE
         IF (film) THEN
            DO 250 k = 1,nvac
               DO 230 j = 1,odi%nq2 - 1
                  DO 220 i = 1,nmzxy
                     vxy(i,j,k,jspins) = vxy(i,j,k,1)
  220             CONTINUE
  230          CONTINUE
               DO 240 i = 1,nmz
                  vz(i,k,jspins) = vz(i,k,1)
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE
         CLOSE (nt)
      END IF
      IF (total) THEN
         OPEN (11,file='potcoul',form='unformatted',status='unknown')
         DO js = 1,jspins
           DO i = 1,nq3                          ! to enable a GW calculation,
             vpw_w(i,js)=vpw_w(i,js)/nstr(i)     ! the PW-coulomb part is not
           ENDDO                                 ! used otherwise anyway.
         ENDDO
         CALL wrtdop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >               jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >               invs,invs2,film,delz,z1,dx,rmt,zatom,
     >               nlh,jri,ntypsd,ntypsy,namat,11,
     >               iop,dop,iter,vr,vpw_w,vz,vxy,name)
         DO js = 1,jspins
           DO i = 1,nq3
             vpw_w(i,js)=vpw_w(i,js)*nstr(i)
           ENDDO
         ENDDO
         CLOSE(11)
      END IF
      IF (plpot) THEN
        OPEN (11,file='potcoul_pl',form='unformatted',status='unknown')
        CALL wrtdop(
     >              jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >              jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >              invs,invs2,film,delz,z1,dx,rmt,zatom,
     >              nlh,jri,ntypsd,ntypsy,namat,11,
     >              iop,dop,iter,vr,vpw,vz,vxy,name)
        CLOSE(11)
      ENDIF

c     ******** exchange correlation potential******************
c+ta
c     rhmn: rho-min.
c     ichsmrg: i-charge-small-region.
c          0(not watched)
c          1(in muffin-tin), 2(interstitial), 3(warped-vac),4(vacuum)

      ichsmrg=0
      rhmn=1.e+10
c-ta
      excpw(:) = cmplx(0.,0.)
      excz(:,:) = 0.0
      excxy(:,:,:) = cmplx(0.,0.)
      excr(:,:,:) = 0.0

c     ---> vacuum region
      IF (film) THEN
         CALL cpu_time(time1)

         ifftd2 = 9*k1d*k2d
         if (odi%d1) ifftd2 = 9*k3d*odi%M

         IF ((igrd == 0).AND.(icorr /= -1)) THEN  ! LDA

          IF (.not.odi%d1) THEN

           CALL vvacxc(
     >                 k1d,k2d,nmzd,nmzxyd,n2d,jspd,ifftd2,
     >                 icorr,total,krla,nmzxy,jspins,nq2,nmz,nstr2,
     >                 rhtxy,rht,cdomvxy,cdomvz,l_noco,
     >                 kimax2,igfft2,pgfft2,nvac,
     X                 vxy,vz,
     <                 excxy,excz)

          ELSE

           CALL vvacxc(
     >                 k3d,odi%M,nmzd,nmzxyd,odi%n2d,jspd,ifftd2,
     >                 icorr,total,krla,nmzxy,jspins,odi%nq2,nmz,
     >                 odi%nst2,rhtxy,rht,cdomvxy,cdomvz,l_noco,
     >                 odi%kimax2,odl%igf,odl%pgf,nvac,
     X                 vxy,vz,
     <                 excxy,excz)

          ENDIF

         ELSE      ! GGA

          IF (odi%d1) THEN

           CALL vvacxcg(
     >          k3d,odi%M,nmzd,nmzxyd,odi%n2d,jspd,odi%nn2d,
     >          ifftd2,nvac,z1,delz,icorr,total,nmzxy,jspins,
     >          odi%nq2,bmat,nmz,odi%nst2,
     >          odg%pgfx,odg%pgfy,odg%pgfxx,odg%pgfyy,odg%pgfxy,
     >          igrd,ndvgrd,isprsv,idsprsv,chng,sprsv,ichsmrg,
     >          rhtxy,rht,cdomvxy,cdomvz,l_noco,l_ss,qss,
     >          odi%kimax2,odl%igf,odl%pgf,odi,
     X          vxy,vz,idsprs,rhmn,
     <          excxy,excz)

          ELSE

           CALL vvacxcg(
     >          k1d,k2d,nmzd,nmzxyd,n2d,jspd,nn2d,ifftd2,nvac,z1,
     >          delz,icorr,total,nmzxy,jspins,nq2,bmat,nmz,nstr2,
     >          pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy,
     >          igrd,ndvgrd,isprsv,idsprsv,chng,sprsv,ichsmrg,
     >          rhtxy,rht,cdomvxy,cdomvz,l_noco,l_ss,qss,
     >          kimax2,igfft2,pgfft2,odi,
     X          vxy,vz,idsprs,rhmn,
     <          excxy,excz)

          END IF
            
         END IF

         CALL cpu_time(time2)
         CALL outtime('generation of ex.-cor.potential in vacuum:',
     +                       time2-time1)
c+odim
      END IF
c     ----------------------------------------
c     ---> interstitial region
      CALL cpu_time(time1)

      ifftd=27*k1d*k2d*k3d

      IF ( (.not. lwb) .or. ( (igrd == 0) .and. (icorr /= -1) ) ) THEN
      ! no White-Bird-trick

        ifftxc3d = kxc1d*kxc2d*kxc3d

        IF ( (igrd == 0) .and. (icorr /= -1) ) THEN
        ! LDA

          CALL visxc(
     >            ifftd,k1d,k2d,k3d,n3d,jspd,
     >            qpw,cdom,l_noco,
     >            kimax,igfft,pgfft,ufft,
     >            icorr,total,krla,jspins,nq3,nstr,
     X            vpw,vpw_w,
     <            excpw)

        ELSE ! GGA

          CALL visxcg(
     >                ifftd,k1d,k2d,k3d,n3d,jspd,nop,
     >                kxc1d,kxc2d,kxc3d,ifftxc3d,
     >                kv3,nop2,mrot,bmat,kxc1_fft,kxc2_fft,
     >                kxc3_fft,nxc3_fft,kmxxc_fft,
     >                qpw,cdom,kimax,igfft,pgfft,ufft,
     >                icorr,total,jspins,nq3,nstr,
     >                igrd,ndvgrd,idsprs,isprsv,l_noco,l_ss,qss,
     >                idsprsi,chng,sprsv,lwb,rhmn,ichsmrg,
     X                vpw,vpw_w,
     <                excpw)

        END IF

        DEALLOCATE ( cdom,cdomvz,cdomvxy )

      ELSE
      ! White-Bird-trick

          WRITE(6,'('' W+B trick cancelled out. visxcwb uses at'',
     +              '' present common block cpgft3. visxcwb needs'',/,
     +              '' to be reprogrammed according to visxcg.f'')')
          STOP 'visxcwb'
csb       CALL visxcwb(
csb  >                 qpw,kimax,igfft,pgfft,ufft,
csb  >                 icorr,total,krla,
csb  >                 igrd,ndvgrd,idsprs,isprsv,
csb  >                 idsprsi,chng,sprsv,lwb,rhmn,ichsmrg,
csb  =                 vpw,vpw_w,
csb  <                 excpw)

      END IF
c
c --> on output vpw_w contains the warped effective potential and
c               excpw the warped XC-energy density
c
      CALL cpu_time(time2)
      CALL outtime('generation of ex.-cor.potential in interstitial:',
     +                       time2-time1)
c+gu
c           ---> evaluate the interstitial average potential
       DO 360 js = 1,jspins
         vbar(js) = vpw_w(1,js)*omtil/volint
  360 CONTINUE
c-gu
      WRITE (6,FMT=8040) (vbar(js),js=1,jspins)
      WRITE (16,FMT=8040) (vbar(js),js=1,jspins)
 8040 FORMAT (/,5x,'interstitial potential average (vbar) =',2f10.6)
c
c     ------------------------------------------
c     ----> muffin tin spheres region
      CALL cpu_time(time1)

      IF ((igrd.EQ.0).AND.(icorr.NE.-1)) THEN
         CALL vmtxc(
     >              jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,nspd,
     >              clnu,mlh,nmem,llh,nlh,rmsh,ntypsy,jri,natd,neq,
     >              rho,icorr,total,krla,ntype,jspins,nsymt,
     X              vr,
     <              excr)
      ELSEIF ((igrd.GT.0).OR.(icorr.EQ.-1)) THEN
         CALL vmtxcg(
     >               jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,nspd,
     >               clnu,mlh,nmem,llh,nlh,rmsh,ntypsy,jri,dx,natd,neq,
     >               rho,icorr,total,krla,ntype,jspins,nsymt,
     >               igrd,ndvgrd,isprsv,idsprs0,idsprsl,chng,sprsv,
     X               vr,idsprs,rhmn,ichsmrg,
     <               excr)
      ELSE
        STOP 'something wrong with igrd before vmtxc'
      ENDIF

      CALL cpu_time(time2)
      CALL outtime('generation of ex.-cor.potential in mt:',time2-time1)
c     ------------------------------------------
c     ---> check continuity of total potential

      IF (vchk) THEN
c           ----> vacuum boundaries
         IF (film .AND. .NOT.odi%d1) THEN
            npd = min(nspd,25)
            CALL points(xp,npd)
            DO ivac = 1,nvac
              signum = 3. - 2.*ivac
              DO j = 1,npd
                 xp(3,j) = signum*z1/amat(3,3)
              ENDDO
              CALL checkdop(
     >                    xp,npd,0,0,ivac,1,1,.false.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    vpw,vr,vxy,vz,odi,ods)
            ENDDO ! ivac = 1,nvac
         ELSEIF (odi%d1) THEN
c-odim
            npd = min(nspd,25)
            CALL cylpts(xp,npd,z1)
c           DO j = 1,npd
c              xp(1,j) = xp(1,j)/amat(1,1)
c              xp(2,j) = xp(2,j)/amat(2,2)
c           ENDDO
            CALL checkdop(
     >           xp,npd,0,0,nvac,1,1,.false.,nspd,jmtd,
     >           memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >           jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >           nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >           tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >           vpw,vr,vxy,vz,odi,ods)
c+odim
         END IF
c           ----> m.t. boundaries
         nat = 1
         DO n = 1, ntype
            CALL sphpts(xp,nspd,rmt(n),pos(1,nat))
            CALL checkdop(
     >                    xp,nspd,n,nat,0,-1,1,.false.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    vpw,vr,vxy,vz,odi,ods)
            nat = nat + neq(n)
         ENDDO ! n = 1, ntype
      END IF

!pg      vlm(:,:,:) = cmplx(0.0,0.0)
!pg      nat = 1
!pg      DO n = 1, ntype
!pg         ns = ntypsy(nat)
!pg         jm = jri(n)
!pg         DO nl = 0, nlh(ns)
!pg            l = llh(nl,ns)
!pg            DO j = 1, jm
!pg               f(j) = (rmsh(j,n)**l)*vr(j,nl,n,1)
!pg            ENDDO
!pg            CALL intgr3(f,rmsh(1,n),dx(n),jm,fint)
!pg            DO mb = 1, nmem(nl,ns)
!pg               m = mlh(mb,nl,ns)
!pg               vlm(m,l,n) = vlm(m,l,n) + clnu(mb,nl,ns)*fint
!pg            ENDDO
!pg         ENDDO
!pg         nat = nat + neq(n)
!pg      ENDDO
!pg      fint = sqrt(2*pimach())/0.529177  ! conversion factor
!pg      write(6,*) 'atom 1'
!pg      vlm(:,1,1) = vlm(:,1,1) * fint
!pg      write(6,'(a6,2f12.6,a4)') 'E_x = ',vlm(-1,1,1)-vlm(1,1,1),' V/A'
!pg      write(6,'(a6,2f12.6,a4)') 'E_y = ',vlm(-1,1,1)+vlm(1,1,1),' V/A'
!pg      write(6,'(a6,2f12.6,a4)') 'E_z = ',vlm(0,1,1),' V/A'
!pg      write(6,*) 'atom 2'
!pg      vlm(:,1,2) = vlm(:,1,2) * fint
!pg      write(6,'(a6,2f12.6,a4)') 'E_x = ',vlm(-1,1,2)-vlm(1,1,2),' V/A'
!pg      write(6,'(a6,2f12.6,a4)') 'E_y = ',vlm(-1,1,2)+vlm(1,1,2),' V/A'
!pg      write(6,'(a6,2f12.6,a4)') 'E_z = ',vlm(0,1,2),' V/A'
!pg      STOP
       
      CALL pot_mod(
     >             jmtd,nlhd,ntypd,nmzd,nmzxyd,n2d,n3d,natd,
     >             ntypsd,jspins,film,nq2,nq3,ntype,
     >             neq,nvac,nmz,nmzxy,jri,ntypsy,nlh,
     X             vr,vxy,vz,vpw,vpw_w)
C
C============TOTAL======================================
C
      IF (total) THEN

        IF (l_noco) THEN ! load qpw,rht,rhtxy from 'cdn'-file
          nt = 70
          OPEN (nt,file='cdn',form='unformatted',status='old')
          CALL loddop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >               jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >               nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <               iop,dop,iter,rho,qpw,rht,rhtxy,name)
          CLOSE (nt)
        ENDIF
C
C     CALCULATE THE INTEGRAL OF n1*Veff1 + n2*Veff2
C     Veff = Vcoulomb + Vxc
C
        te_veff = 0.0
        DO 370 js = 1,jspins
           WRITE (6,FMT=8050) js
           WRITE (16,FMT=8050) js
 8050 FORMAT (/,10x,'density-effective potential integrals for spin ',
     +          i2,/)

           CALL int_nv(
     >             n3d,n2d,nmzxyd,nmzd,jmtd,ntypd,nlhd,ntypsd,
     >             qpw(1,js),vpw_w(1,js),
     >             rhtxy(1,1,1,js),vxy(1,1,1,js),
     >             rht(1,1,js),vz(1,1,js),
     >             rho(1,0,1,js),vr(1,0,1,js),
     >             area,delz,omtil,rmsh,dx,
     >             jri,neq,nlh,nstr2,ntypsy,natd,
     >             nq3,
     >             nmz,nmzxy,nq2,ntype,nvac,
     >             zrfs,film,invs,z1,odi,
     <             te_veff)
 370     CONTINUE
      WRITE (6,FMT=8060) te_veff
      WRITE (16,FMT=8060) te_veff
 8060 FORMAT (/,10x,'total density-effective potential integral :',
     +     t40,f20.10)
C
C     CALCULATE THE INTEGRAL OF n*exc
C   
c     ---> perform spin summation of charge densities
c     ---> for the calculation of Exc
      IF (jspins.EQ.2) THEN
         nat = 1
         DO 400 n = 1,ntype
            DO 390 lh = 0,nlh(ntypsy(nat))
               DO 380 i = 1,jri(n)
                  rho(i,lh,n,1) = rho(i,lh,n,1) + rho(i,lh,n,jspins)
  380          CONTINUE
  390       CONTINUE
            nat = nat + neq(n)
  400    CONTINUE
         DO 410 i = 1,nq3
            qpw(i,1) = qpw(i,1) + qpw(i,jspins)
  410    CONTINUE
         IF (film) THEN
            DO 450 ivac = 1,nvac
               DO 430 j = 1,odi%nq2 - 1
                  DO 420 i = 1,nmzxy
                     rhtxy(i,j,ivac,1) = rhtxy(i,j,ivac,1) +
     +                                   rhtxy(i,j,ivac,jspins)
  420             CONTINUE
  430          CONTINUE
               DO 440 i = 1,nmz
                  rht(i,ivac,1) = rht(i,ivac,1) + rht(i,ivac,jspins)
  440          CONTINUE
  450       CONTINUE
         END IF
      END IF
      WRITE (6,FMT=8070)
      WRITE (16,FMT=8070)
 8070 FORMAT (/,10x,'charge density-energy density integrals',/)

      te_exc = 0.0
      CALL int_nv(
     >             n3d,n2d,nmzxyd,nmzd,jmtd,ntypd,nlhd,ntypsd,
     >             qpw(1,1),excpw(1),
     >             rhtxy,excxy,
     >             rht,excz,
     >             rho,excr,
     >             area,delz,omtil,rmsh,dx,
     >             jri,neq,nlh,nstr2,ntypsy,natd,
     >             nq3,
     >             nmz,nmzxy,nq2,ntype,nvac,
     >             zrfs,film,invs,z1,odi,
     <             te_exc)
      WRITE (6,FMT=8080) te_exc
      WRITE (16,FMT=8080) te_exc
      
 8080 FORMAT (/,10x,'total charge density-energy density integral :',
     +        t40,f20.10)

      END IF
      IF (plpot) THEN
        OPEN (11,file='pottot_pl',form='unformatted',status='unknown')
        CALL wrtdop(
     >              jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >              jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >              invs,invs2,film,delz,z1,dx,rmt,zatom,
     >              nlh,jri,ntypsd,ntypsy,namat,11,
     >              iop,dop,iter,vr,vpw,vz,vxy,name)
        CLOSE(11)
      ENDIF
C
C==========END TOTAL============================================
C
c           ---> store v(l=0) component as r*v(l=0)/sqrt(4pi)

      DO js = 1,jspins

         l_pottot = .FALSE.                  ! adds a B-field, if file
         IF (jspd.eq.2) THEN                       !   mfee is present
           INQUIRE (file='mfee',exist=l_pottot)
           IF (l_pottot) THEN
             OPEN (88,file='mfee',form='formatted',status='unknown')
             REWIND 88
             DO n=1,ntype
               READ (88,'(i2,1x,f8.5)') i,mfie
               WRITE (*,*) 'type,field:',i,mfie
               IF (i.NE.n) STOP 'wrong types in mfee'
                 DO j = 1,jri(n)
                   IF (js.eq.1) then
                    vr(j,0,n,js) = vr(j,0,n,js) - mfie/2.
                   ELSE
                    vr(j,0,n,js) = vr(j,0,n,js) + mfie/2.
                   ENDIF 
                 ENDDO
             ENDDO
             CLOSE (88)
           ENDIF
         ENDIf


         DO n = 1,ntype
            DO j = 1,jri(n)
              vr(j,0,n,js) = rmsh(j,n)*vr(j,0,n,js)/sfp
            ENDDO 
         ENDDO 

      ENDDO     ! js =1,jspins

#ifdef CPP_GF
      OPEN(12,file='gf_intpot',form='unformatted',status='replace')
          CALL wrtdop(
     >               jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >               jspins,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     >               invs,invs2,film,delz,z1,dx,rmt,zatom,
     >               nlh,jri,ntypsd,ntypsy,namat,12,
     >               'potentia','l       ',iter,vr,vpw,vz,vxy,name)
c      write(12) n3d
c      write (12) vpw
       CLOSE(12)
#endif

      IF ((.not.reap).OR.(l_noco)) THEN
         IF (total) THEN
           OPEN (9,file='nrp',form='unformatted',status='unknown')
           REWIND 9
         ENDIF
         CALL wrtdop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >               jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >               invs,invs2,film,delz,z1,dx,rmt,zatom,
     >               nlh,jri,ntypsd,ntypsy,namat,9,
     >               'potentia','l       ',iter,vr,vpw,vz,vxy,name)
         CLOSE(9)
      ENDIF

c     **************** reanalyze vpw *************************
c                        call cpu_time(cp0)
      IF (total) THEN
c     ----->write potential to file 8
c -> the following procedure is required in order to run
c    correctly on the helga parallel cluster in Hamburg end 2005
c                                                Paolo & YM
         l_pottot = .FALSE.
         INQUIRE (file='pottot',exist=l_pottot)
         IF (l_pottot) THEN
            OPEN (8,file='pottot',form='unformatted',status='unknown')
            CLOSE (8,status='delete')
            WRITE(6,*) 'vgen: pottot deleted'
         ENDIF
         OPEN (8,file='pottot',form='unformatted',status='unknown')
         REWIND 8
         DO js=1,jspins
           DO i=1,nq3
             vpw_w(i,js)=vpw_w(i,js)/nstr(i)
           ENDDO
         ENDDO
         CALL wrtdop(
     >               jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >               jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >               invs,invs2,film,delz,z1,dx,rmt,zatom,
     >               nlh,jri,ntypsd,ntypsy,namat,8,
     >               'potentia','l       ',iter,vr,vpw_w,vz,vxy,name) ! vpw_w
        CLOSE(8)
      END IF

      ENDIF ! irank == 0

      DEALLOCATE (alphm,excpw,excxy,vbar,af1,bf1,xp,rho,rht,
     +            qpw,rhtxy,vr,vz,vpw,vxy,excz,excr,vpw_w,psq)

      END SUBROUTINE vgen
      END MODULE m_vgen
