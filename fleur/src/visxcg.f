      MODULE m_visxcg
c     ******************************************************
c     subroutine generates the exchange-correlation potential
c     in the interstitial region    c.l.fu
c     including gradient corrections. t.a. 1996.
c     ******************************************************
      CONTAINS
      SUBROUTINE visxcg(
     >                  ifftd,k1d,k2d,k3d,n3d,jspd,nop,
     >                  kxc1d,kxc2d,kxc3d,ifftxc3d,
     >                  kv3,nop2,mrot,bmat,kxc1_fft,kxc2_fft,
     >                  kxc3_fft,nxc3_fft,kmxxc_fft,
     >                  qpw,cdom,kimax,igfft,pgfft,ufft,
     >                  icorr,total,jspins,ng3,nstr,
     >                  igrd,ndvgrd,idsprs,isprsv,l_noco,l_ss,qss,
     >                  idsprsi,chng,sprsv,lwb,rhmn,ichsmrg,
     X                  vpw,vpw_w,
     <                  excpw)

c     ******************************************************
c     instead of visxcor.f: the different exchange-correlation
c     potentials defined through the key icorr are called through
c     the driver subroutine vxcallg.f,for the energy density - excallg
c     subroutines vectorized
c     ** r.pentcheva 22.01.96
c     *********************************************************
c     in case of total = .true. calculates the ex-corr. energy
c     density
c     ** r.pentcheva 08.05.96
c     ******************************************************************

      USE m_grdrsis
      USE m_prpxcfftmap
      USE m_mkgxyz3
      USE m_xcallg, ONLY : vxcallg,excallg
      USE m_fft3d
      USE m_fft3dxc
      USE m_set
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ifftd,k1d,k2d,k3d,n3d,jspd,nop
      INTEGER, INTENT (IN) :: kxc1d,kxc2d,kxc3d,ifftxc3d
      INTEGER, INTENT (IN) :: icorr,jspins,ng3,igrd,ndvgrd 
      LOGICAL, INTENT (IN) :: total,l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nstr(n3d)

      REAL rhmn,rhmni,chng,d_15,sprsv
      INTEGER ichsmrg,idsprs,isprsv,idsprsi
c     lwb: if true, white-bird trick.
c     lwbc: l-white-bird-current.
      LOGICAL lwb,lwbc
c
c----->  fft  information  for xc potential + energy
c
      INTEGER, INTENT (IN) :: nxc3_fft,kmxxc_fft
      INTEGER, INTENT (IN) :: kxc1_fft,kxc2_fft,kxc3_fft
      INTEGER, ALLOCATABLE :: igxc_fft(:)
      REAL,    ALLOCATABLE :: gxc_fft(:,:)
c
c----->  fft  information  for convolution of step-function * potential
c
      INTEGER, INTENT (IN) :: kimax
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
c
c-----> stars and symmetry information
c
      INTEGER, INTENT (IN) :: kv3(3,n3d),nop2,mrot(3,3,nop)
      REAL,    INTENT (IN) :: bmat(3,3)
c
c-----> charge density, potential and energy density
c
      COMPLEX, INTENT (IN) :: qpw(n3d,jspd)
      COMPLEX, INTENT (OUT) :: excpw(n3d)
      COMPLEX, INTENT (INOUT) ::vpw(n3d,jspd),vpw_w(n3d,jspd),cdom(n3d)
c     ..
c     .. Local Scalars ..
      INTEGER :: i,ig,k,js,nt,ivec(n3d),idsprssv,ifftxc3,idm,jdm,ndm
      COMPLEX :: ci
      REAL    :: rhotot
c     ..
c     .. Local Arrays ..
      COMPLEX :: fg3(n3d),cqpw(n3d,jspins),ccdom(n3d)
      REAL, ALLOCATABLE :: ph_wrk(:),bf3(:)
      REAL, ALLOCATABLE :: rho(:,:),rhd1(:,:,:),rhd2(:,:,:)
      REAL, ALLOCATABLE :: mx(:),my(:)
      REAL, ALLOCATABLE :: magmom(:),dmagmom(:,:),ddmagmom(:,:,:) 
c 
      REAL, ALLOCATABLE :: vxc(:,:),exc(:),vcon(:) 
      REAL, ALLOCATABLE :: agr(:),agru(:),agrd(:)
      REAL, ALLOCATABLE :: g2r(:),g2ru(:),g2rd(:)
      REAL, ALLOCATABLE :: gggr(:),gggru(:),gggrd(:)
      REAL, ALLOCATABLE :: gzgr(:)

c     .. unused input (needed for other noco GGA-implementations) ..
      LOGICAL, INTENT (IN) :: l_ss    
      REAL,    INTENT (IN) :: qss(3) 

cta+
c.....------------------------------------------------------------------
c------->          abbreviations
c
c     ph_wrk: work array containing phase * g_x,gy...... 
c     qpw   : charge density stored as stars
c     rho   : charge density stored in real space
c     vxc   : exchange-correlation potential in real space
c     exc   : exchange-correlation energy density in real space
c     kxc1d  : dimension of the charge density fft box in the pos. domain
c     kxc2d  : defined in dimens.f program (subroutine apws). 1,2,3 indic
c     kxc3d  ; a_1, a_2, a_3 directions.
c     kq(i) : i=1,2,3 actual length of the fft-box for which fft is done
c     nstr  : number of members (arms) of reciprocal lattice (g) vector
c             of each star.
c     nxc3_fft: number of stars in the  charge density  fft-box
c     ng3   : number of 3 dim. stars in the charge density sphere define
c             by gmax
c     kmxxc_fft: number of g-vectors forming the nxc3_fft stars in the
c               charge density or xc-density sphere
c     kimax : number of g-vectors forming the ng3 stars in the gmax-sphe
c     ifftxc3d: elements (g-vectors) in the charge density  fft-box
c     igfft : pointer from the g-sphere (stored as stars) to fft-grid
c             and     from fft-grid to g-sphere (stored as stars)
c     pgfft : contains the phases of the g-vectors of sph.
c     isn   : isn = +1, fft transform for g-space to r-space
c             isn = -1, vice versa
c
c-------------------------------------------------------------------
c
c---> set up pointer for backtransformation of from g-vector in
c     positive domain of xc density fftbox into stars.
c     also the x,y,z components of the g-vectors are set up to calculate
c     derivatives.
c     in principle this can also be done in main program once.
c     it is done here to save memory.
c

      ALLOCATE ( igxc_fft(0:ifftxc3d-1),gxc_fft(0:ifftxc3d-1,3) )
      CALL prp_xcfft_map(
     >                   n3d,kxc1d,kxc2d,kxc3d,nop,
     >                   kv3,nop2,mrot,bmat,
     >                   kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,
     <                   igxc_fft,gxc_fft)
c
      ifftxc3=kxc1_fft*kxc2_fft*kxc3_fft
      idsprssv=idsprs
      idsprs=idsprsi

      lwbc=lwb

      IF (ng3.GT.n3d) THEN
        WRITE(6,'(/'' ng3.gt.n3d. ng3,n3d='',2i6)') ng3,n3d
        STOP 'visxcg: ng3.gt.n3d'
      ENDIF

      d_15=1.e-15
c
      ci=cmplx(0.,1.)
c
c Allocate arrays
c ff
      ALLOCATE( bf3(0:ifftd-1),ph_wrk(0:ifftxc3d-1),  
     +          rho(0:ifftxc3d-1,jspd),rhd1(0:ifftxc3d-1,jspd,3),
     +          rhd2(0:ifftxc3d-1,jspd,6) )
      IF (l_noco)  THEN
        ALLOCATE( mx(0:ifftxc3-1),my(0:ifftxc3-1),
     +            magmom(0:ifftxc3-1),  
     +            dmagmom(0:ifftxc3-1,3),ddmagmom(0:ifftxc3-1,3,3) )
      END IF  


c-->     transform charge density to real space

      DO js=1,jspins
        CALL fft3dxc(
     <              rho(0,js),bf3,
     >              qpw(1,js),
     >              kxc1_fft,kxc2_fft,kxc3_fft,
     >              nxc3_fft,kmxxc_fft,+1,
     >              igfft(0,1),igxc_fft,pgfft,nstr)
      END DO

      IF (l_noco) THEN  

c       for off-diagonal parts the same
        CALL fft3dxc(
     <               mx,my,
     >               cdom,
     >               kxc1_fft,kxc2_fft,kxc3_fft,
     >               nxc3_fft,kmxxc_fft,+1,
     >               igfft(0,1),igxc_fft,pgfft,nstr)

        DO i=0,ifftxc3-1 
          rhotot= 0.5*( rho(i,1) + rho(i,2) )
          magmom(i)= SQRT(  (0.5*(rho(i,1)-rho(i,2)))**2 
     +                    + mx(i)**2 + my(i)**2 )
          rho(i,1)= rhotot+magmom(i)
          rho(i,2)= rhotot-magmom(i)
        END DO 

      ENDIF

      IF (igrd.EQ.0) GOTO 100  

! In collinear calculations all derivatives are calculated in g-spce,
! in non-collinear calculations the derivatives of |m| are calculated in real space. 

c-->   for d(rho)/d(x,y,z) = rhd1(:,:,idm) (idm=1,2,3).
c
c         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z).

      DO js= 1,jspins
        DO i = 1,ng3
          cqpw(i,js)= ci*qpw(i,js)
        END DO
      END DO

      DO idm=1,3

        DO ig = 0 , kmxxc_fft - 1
          ph_wrk(ig) = pgfft(ig) * gxc_fft(ig,idm)
        END DO

        DO js=1,jspins
          CALL fft3dxc(
     <           rhd1(0,js,idm),bf3,
     >           cqpw(1,js),
     >           kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,+1,
     >           igfft(0,1),igxc_fft,ph_wrk,nstr)
        END DO

      END DO

      IF (l_noco) THEN

        CALL grdrsis(
     >           magmom,bmat,kxc1_fft,kxc2_fft,kxc3_fft,ndvgrd,
     <           dmagmom )

        DO i=0,ifftxc3-1
          DO idm=1,3
            rhotot= rhd1(i,1,idm)/2.+rhd1(i,2,idm)/2.
            rhd1(i,1,idm)= rhotot+dmagmom(i,idm) 
            rhd1(i,2,idm)= rhotot-dmagmom(i,idm) 
          END DO  
        END DO 

      END IF  

      IF (lwbc) GOTO 100 

c-->   for dd(rho)/d(xx,xy,yy,zx,yz,zz) = rhd2(:,:,idm) (idm=1,2,3,4,5,6)
c
c         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z) * g_(x,y,z)

      DO i = 1,ng3
        DO js=1,jspins 
          cqpw(i,js)= -qpw(i,js)
        END DO
      END DO

      ndm = 0
      DO idm = 1,3
        DO jdm = 1,idm
          ndm = ndm + 1

          DO ig = 0 , kmxxc_fft-1
            ph_wrk(ig) = pgfft(ig)*gxc_fft(ig,idm)*gxc_fft(ig,jdm)
          ENDDO

          DO js=1,jspins
            CALL fft3dxc(
     <                rhd2(0,js,ndm),bf3,
     >                cqpw(1,js),
     >                kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,+1,
     >                igfft(0,1),igxc_fft,ph_wrk,nstr)
          END DO
        END DO ! jdm 
      END DO   ! idm 

      IF (l_noco) THEN

        DO idm = 1,3
          CALL grdrsis(
     >           dmagmom(0,idm),bmat,kxc1_fft,kxc2_fft,kxc3_fft,ndvgrd,
     <           ddmagmom(0,1,idm) )
        END DO 

        ndm= 0
        DO idm = 1,3
          DO jdm = 1,idm
            ndm = ndm + 1  

            DO i=0,ifftxc3-1
              rhotot= rhd2(i,1,ndm)/2.+rhd2(i,2,ndm)/2.
              rhd2(i,1,ndm)= rhotot +
     &         ( ddmagmom(i,jdm,idm) + ddmagmom(i,idm,jdm) )/2. 
              rhd2(i,2,ndm)= rhotot -
     &         ( ddmagmom(i,jdm,idm) + ddmagmom(i,idm,jdm) )/2. 
            END DO

          ENDDO !jdm
        ENDDO   !idm 
 
      END IF  

  100 CONTINUE


      DEALLOCATE ( ph_wrk )
      IF (l_noco) THEN 
        DEALLOCATE(mx,my,magmom,dmagmom,ddmagmom) 
      END IF  
c
      DO js=1,jspins 
        DO i=0,ifftxc3-1
          rho(i,js)=max(rho(i,js),d_15)
        ENDDO
      END DO

      CALL dset(ifftd,0.0,bf3,1)
c
c allocate the other arrays 
c
      ALLOCATE (agr(0:ifftxc3d-1),agru(0:ifftxc3d-1),agrd(0:ifftxc3d-1),
     +          g2r(0:ifftxc3d-1),g2ru(0:ifftxc3d-1),g2rd(0:ifftxc3d-1),
     +       gggr(0:ifftxc3d-1),gggru(0:ifftxc3d-1),gggrd(0:ifftxc3d-1),
     +       gzgr(0:ifftxc3d-1))

c
c     calculate the quantities such as abs(grad(rho)),.. used in
c     evaluating the gradient contributions to potential and energy.
c
      CALL mkgxyz3
     >            (igrd,ifftxc3d,jspd,ifftxc3,jspins,rho,
     >             rhd1(0,1,1),rhd1(0,1,2),rhd1(0,1,3),
     >             rhd2(0,1,1),rhd2(0,1,3),rhd2(0,1,6),
     >             rhd2(0,1,5),rhd2(0,1,4),rhd2(0,1,2),
     <             agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,
     <             gzgr)

      DEALLOCATE ( rhd1,rhd2 )
      ALLOCATE ( vxc(0:ifftxc3d-1,jspd) )
c
c     calculate the exchange-correlation potential in  real space
c
       nt=ifftxc3
c
c      rhmni: rho_minimum_interstitial.

       rhmni=10.e+10

       DO js=1,jspins
        DO i=0,ifftxc3-1
          rho(i,js)=max(rho(i,js),d_15)
          rhmni=min(rhmni,rho(i,js))
        ENDDO
       ENDDO

       IF (rhmni.lt.rhmn) THEN
         rhmn=rhmni
         ichsmrg=2
       ENDIF

       IF (rhmn.lt.chng) then
         WRITE(6,'(/'' rhmn.lt.chng in visxc. rhmn,chng='',
     +     2d9.2)') rhmn,chng
c         STOP 'visxcg: rhmn.lt.chng'
       ENDIF

       CALL vxcallg(icorr,lwbc,jspins,ifftxc3d,nt,rho,
     >              agr,agru,agrd,g2r,g2ru,g2rd,
     &              gggr,gggru,gggrd,gzgr,
     <              vxc,
     >              idsprs,isprsv,sprsv)
c
c----> back fft to g space
c----> perform back  fft transform: vxc(r) --> vxc(star)
c
      DO i = 1,ng3
        ivec(i)=1
      ENDDO

      DO js = 1,jspins
        CALL dset(ifftxc3d,0.0,bf3(0),1)
        CALL fft3dxc(
     >               vxc(0,js),bf3,
     <               fg3,
     >               kxc1_fft,kxc2_fft,kxc3_fft,
     >               nxc3_fft,kmxxc_fft,-1,
     >               igfft(0,1),igxc_fft,pgfft,nstr)
c
        DO k = 1,nxc3_fft
           vpw(k,js) = vpw(k,js) + fg3(k)
        ENDDO

c
c====>  INCLUDING TOTAL ENERGY
c
        IF (total) THEN
c
c----> Perform fft transform: vxc(star) --> vxc(r) 
c     !Use large fft mesh for convolution
c
           IF (nxc3_fft.LT.ng3) THEN
              CALL cset((ng3-nxc3_fft),(0.0,0.0),fg3(nxc3_fft+1),1)
           ENDIF
           ALLOCATE ( vcon(0:ifftd-1) )
           CALL fft3d(
     <                vcon(0),bf3,
     >                fg3,
     >                k1d,k2d,k3d,ng3,kimax,+1,
     >                igfft,pgfft,nstr)
c
c----> Convolute with step function
c
           DO i=0,ifftd-1
              vcon(i)=ufft(i)*vcon(i)
           ENDDO

           CALL dset(ifftd,0.0,bf3(0),1)
           CALL fft3d(
     >                vcon(0),bf3,
     <                fg3,
     >                k1d,k2d,k3d,ng3,kimax,-1,
     >                igfft,pgfft,ivec)
           DEALLOCATE ( vcon )
c
c----> add to warped coulomb potential
c
           DO k = 1,ng3
              vpw_w(k,js) = vpw_w(k,js) + fg3(k)
           ENDDO

         ENDIF
      ENDDO
      DEALLOCATE ( vxc )
c
c     calculate the ex.-cor energy density in real space
c
      IF (total) THEN
         ALLOCATE ( exc(0:ifftxc3d-1) )
         CALL excallg(icorr,lwbc,jspins,ifftxc3d,nt,rho,
     >               agr,agru,agrd,g2r,g2ru,g2rd,
     &               gggr,gggru,gggrd,gzgr,
     <               exc,
     >               idsprs,isprsv,sprsv)
c
c---->   perform back  fft transform: exc(r) --> exc(star)
c
         CALL dset(ifftxc3d,0.0,bf3(0),1)
         CALL fft3dxc(
     >                exc,bf3,
     <                fg3,
     >                kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,-1,
     >                igfft(0,1),igxc_fft,pgfft,nstr)
         DEALLOCATE ( exc )
c
c---->   Perform fft transform: exc(star) --> exc(r) 
c        !Use large fft mesh for convolution
c
         IF (nxc3_fft.LT.ng3) THEN
            CALL cset((ng3-nxc3_fft),(0.0,0.0),fg3(nxc3_fft+1),1)
         ENDIF
         ALLOCATE ( vcon(0:ifftd-1) )
         CALL fft3d(vcon,bf3,fg3,
     >             k1d,k2d,k3d,ng3,kimax,+1,
     >             igfft,pgfft,nstr)

         DO i=0,ifftd-1
           vcon(i)=ufft(i)*vcon(i)
         ENDDO
c
c         ---> back fft to g space
c
         CALL dset(ifftd,0.0,bf3(0),1)
         CALL fft3d(vcon,bf3,excpw,
     >             k1d,k2d,k3d,ng3,kimax,-1,
     >             igfft,pgfft,ivec)
         DEALLOCATE ( vcon )
c
      ENDIF

      DEALLOCATE ( bf3,rho,igxc_fft,gxc_fft )
      DEALLOCATE ( agr,agru,agrd,g2r,g2ru,g2rd,
     +             gggr,gggru,gggrd,gzgr )

      idsprs=idsprssv

      END SUBROUTINE visxcg
      END MODULE m_visxcg
