      MODULE m_vvacxcg
c-----------------------------------------------------------------------
c     calculates 2-d star function coefficients of exchange-correlation*
c     potential in the vacuum regions and adds them to the corresponding
c     coeffs of the coulomb potential            c.l.fu, r.podloucky   *
c     for the gradient contribution.   t.a. 1996
c-----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE vvacxcg(
     >           k1d,k2d,nmzd,nmzxyd,n2d,jspd,nn2d,ifftd2,nvac,
     >           z1,delz,icorr,total,nmzxy,jspins,ng2,bmat,nmz,nstr2,
     >           pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy,
     >           igrd,ndvgrd,isprsv,idsprsv,chng,sprsv,ichsmrg,
     >           rhtxy,rht,cdomvxy,cdomvz,l_noco,l_ss,qss,
     >           kimax2,igfft2,pgfft2,odi,
     X           vxy,vz,idsprs,rhmn,
     <           excxy,excz)

c-----------------------------------------------------------------------
c     instead of vvacxcor.f: the different exchange-correlation
c     potentials defined through the key icorr are called through
c     the driver subroutine vxcallg.f, subroutines vectorized
c     in case of total = .true. calculates the ex-corr. energy
c     density through the driver subroutine excallg.f
c     ** r.pentcheva 08.05.96
c-----------------------------------------------------------------------
      USE m_od_types
      use m_constants, only : pimach
      USE m_grdrsvac
      USE m_grdchlh
      USE m_mkgz
      USE m_mkgxyz3
      USE m_od_mkgxyz3
      USE m_od_mkgz
      USE m_setmn
      USE m_fft2d
      USE m_set,    ONLY : dset
      USE m_xcallg, ONLY : vxcallg,excallg

      IMPLICIT NONE
c     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,nmzd,nmzxyd,n2d,jspd,nn2d,ifftd2
      INTEGER, INTENT (IN) :: icorr,kimax2,nmzxy,jspins,ng2,nmz
      INTEGER, INTENT (IN) :: igrd,ndvgrd,isprsv,idsprsv,nvac
      REAL,    INTENT (IN) :: sprsv,chng,delz,z1
      LOGICAL, INTENT (IN) :: l_noco,total
      INTEGER, INTENT (INOUT) :: idsprs,ichsmrg
      REAL,    INTENT (INOUT) :: rhmn
c-odim
      TYPE (od_inp), INTENT (INOUT) :: odi
c+odim
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nstr2(n2d)
      INTEGER, INTENT (IN) :: igfft2(0:nn2d-1,2)
      REAL,    INTENT (IN) :: bmat(3,3) 
      REAL,    INTENT (IN) :: pgfft2(0:nn2d-1),pgft2xy(0:nn2d-1)
      REAL,    INTENT (IN) :: pgft2x(0:nn2d-1),pgft2y(0:nn2d-1)
      REAL,    INTENT (IN) :: pgft2xx(0:nn2d-1),pgft2yy(0:nn2d-1)
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd)
      COMPLEX, INTENT (IN) :: rhtxy(nmzxyd,n2d-1,2,jspd)
      COMPLEX, INTENT (IN) :: cdomvz(nmzd,2)
      COMPLEX, INTENT (IN) :: cdomvxy(nmzxyd,n2d-1,2)
      REAL,    INTENT (OUT) :: excz(nmzd,2)
      COMPLEX, INTENT (OUT) :: excxy(nmzxyd,n2d-1,2)
      REAL,    INTENT (INOUT) :: vz(nmzd,2,jspd)
      COMPLEX, INTENT (INOUT) :: vxy(nmzxyd,n2d-1,2,jspd)
c     ..
c     .. Local Scalars ..
      INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip,idsprssv
      REAL    :: rhti,zro,fgz,rhmnv,d_15,bmat1(3,3),rd,tpi
      COMPLEX :: ci 
      LOGICAL :: lwbc              ! if true, white-bird trick.
c     ..
c     .. Local Arrays ..
      REAL, ALLOCATABLE :: af2(:,:),bf2(:),agr(:),agru(:),agrd(:),g2r(:)
      REAL, ALLOCATABLE :: g2ru(:),g2rd(:),gggr(:),gggru(:),gggrd(:)
      REAL, ALLOCATABLE :: gzgr(:),rhdx(:,:),rhdy(:,:),rhdz(:,:)
      REAL, ALLOCATABLE :: rhdxx(:,:),rhdyy(:,:),rhtdz(:,:),rhtdzz(:,:)
      REAL, ALLOCATABLE :: rhdzz(:,:),rhdyz(:,:),rhdzx(:,:),rhdxy(:,:)
      REAL, ALLOCATABLE :: vxc(:,:),exc(:),vxcz(:,:),rxydzr(:),rxydzi(:)
      REAL, ALLOCATABLE :: rxydzzr(:),rxydzzi(:),rhtxyr(:),rhtxyi(:)
      REAL, ALLOCATABLE :: rhtxc(:,:),rhtz(:,:),dummy(:)
      COMPLEX, ALLOCATABLE :: fgxy(:),rxydz(:,:,:),rxydzz(:,:,:),cqpw(:)
c     ..
c     for the noco-case only 
      REAL :: chdens
      REAL, ALLOCATABLE :: magmom(:,:),
     &                     dxmagmom(:),ddxmagmom(:,:),
     &                     dymagmom(:),ddymagmom(:,:), 
     &                     dzmagmom(:,:),ddzmagmom(:,:)
      REAL, ALLOCATABLE :: mx(:),my(:) 
      REAL, ALLOCATABLE :: af2noco(:,:,:) 

c     .. unused input (needed for other noco GGA-implementations) ..
      LOGICAL, INTENT (IN) :: l_ss
      REAL,    INTENT (IN) :: qss(3)

      idsprssv = idsprs
      idsprs   = idsprsv
      lwbc     = .false.
      d_15     = 1.e-15
      zro      = 0.0
      ci       = cmplx(0.,1.)
      nt       = ifftd2
      tpi = 2* pimach()
      if(odi%d1)then
         bmat1(:,:) = 0.
         bmat1(1,1) = bmat(3,3)
         bmat1(2,2) = 1.
      else
         bmat1(:,:) = bmat(:,:)
      endif

      WRITE (6,'(/'' ifftd2,nmz='',2i7)') ifftd2,nmz
      WRITE(6,'('' 9990nmzxy='',2i5)') nmzxy

      ALLOCATE ( rxydz(nmzxyd,n2d-1,jspd),rxydzz(nmzxyd,n2d-1,jspd) )
      ALLOCATE ( rhtdz(nmzd,jspd),rhtdzz(nmzd,jspd) )

      DO ivac=1,nvac 

! the charge density in vacuum is expanded in 2-dim stars on a mesh
! in z-direction. the g||.ne.zero-components expand from 1 to nmzxy
! the g||.eq.zero-components expand from 1 to nmz
! first we calculate vxc in the warping region 

        IF (l_noco) THEN

          ALLOCATE ( magmom(0:ifftd2-1,nmzxy) ) 
          ALLOCATE ( dzmagmom(0:ifftd2-1,nmzxy) ) 
          ALLOCATE ( ddzmagmom(0:ifftd2-1,nmzxy) ) 
          ALLOCATE ( mx(0:ifftd2-1),my(0:ifftd2-1) )
          ALLOCATE ( dxmagmom(0:ifftd2-1),dymagmom(0:ifftd2-1) )
          ALLOCATE ( ddxmagmom(0:ifftd2-1,2),ddymagmom(0:ifftd2-1,2) ) 
          ALLOCATE ( af2noco(0:ifftd2-1,nmzxy,jspd),bf2(0:ifftd2-1) ) 

          ! Transform charge and magnetization to real-space.
          ! In the collinear case that is done later within
          ! another loop over the vacuum-layers in order to 
          ! save memory.

          DO ip=1,nmzxy

            DO js=1,jspins
              CALL fft2d(
     >                   k1d,k2d,n2d,
     <                   af2noco(0,ip,js),bf2,
     >                   rht(ip,ivac,js),0.,rhtxy(ip,1,ivac,js),
     >                   nmzxyd,ng2,kimax2,+1,
     >                   igfft2,pgfft2,nstr2)
            END DO
            CALL fft2d(
     >                 k1d,k2d,n2d,
     <                 mx,my,
     >                 REAL(cdomvz(ip,ivac)),AIMAG(cdomvz(ip,ivac)),
     >                 cdomvxy(ip,1,ivac),
     >                 nmzxyd,ng2,kimax2,+1,
     >                 igfft2,pgfft2,nstr2)

            DO i=0,9*k1d*k2d-1
              magmom(i,ip)= mx(i)**2 + my(i)**2 +
     +          ((af2noco(i,ip,1)-af2noco(i,ip,2))/2.)**2
              magmom(i,ip)= SQRT(magmom(i,ip))
              chdens= af2noco(i,ip,1)/2.+af2noco(i,ip,2)/2.
              af2noco(i,ip,1)= chdens + magmom(i,ip)
              af2noco(i,ip,2)= chdens - magmom(i,ip)
            END DO

          END DO ! ip=1,nmzxy 
          DEALLOCATE ( bf2 )
        END IF ! l_noco 

!      ENDDO    ! ivac
!      DO ivac = 1,nvac

        IF (igrd.GT.0) THEN
          DO js=1,jspins
c
c calculate first (rhtdz) & second (rhtdzz) derivative of rht(1:nmz)
c
            ALLOCATE ( dummy(nmz) )
            CALL grdchlh(
     >                   0,1,nmz,delz,dummy,rht(1,ivac,js),ndvgrd,
     <                   rhtdz(1,js),rhtdzz(1,js))
            DEALLOCATE ( dummy )
            ALLOCATE ( rhtxyr(nmzxy), rhtxyi(nmzxy),dummy(nmzxy) )
            ALLOCATE ( rxydzr(nmzxy), rxydzi(nmzxy) )
            ALLOCATE ( rxydzzr(nmzxy),rxydzzi(nmzxy) )

            DO iq = 1, ng2-1
c
c calculate first (rxydz) & second (rxydzz) derivative of rhtxy:
c
              DO ip=1,nmzxy
                rhtxyr(ip)=rhtxy(ip,iq,ivac,js)
              ENDDO
              CALL grdchlh(
     >                     0,1,nmzxy,delz,dummy,rhtxyr,ndvgrd,
     <                     rxydzr,rxydzzr)

              DO ip=1,nmzxy
                rhtxyi(ip)=aimag(rhtxy(ip,iq,ivac,js))
              ENDDO
              CALL grdchlh(
     >                     0,1,nmzxy,delz,dummy,rhtxyi,ndvgrd,
     <                     rxydzi,rxydzzi)

              DO ip=1,nmzxy
                rxydz(ip,iq,js)=cmplx(rxydzr(ip),rxydzi(ip))
                rxydzz(ip,iq,js)=cmplx(rxydzzr(ip),rxydzzi(ip))
              ENDDO

            ENDDO ! loop over 2D stars (iq)

            DEALLOCATE ( rhtxyr,rhtxyi,rxydzr,rxydzi,rxydzzr,rxydzzi )
            DEALLOCATE ( dummy )

          ENDDO ! jspins

          IF (l_noco) THEN 
!  calculate  dzmagmom = d magmom / d z  and ddzmagmom= d dmagmom / d z 

            ALLOCATE ( rhtxyr(nmzxy),dummy(nmzxy)   )
            ALLOCATE ( rxydzr(nmzxy),rxydzzr(nmzxy) )
            DO i=0,9*k1d*k2d-1 
              DO ip=1,nmzxy
                rhtxyr(ip)=magmom(i,ip)
              ENDDO
              CALL grdchlh(
     >                     0,1,nmzxy,delz,dummy,rhtxyr,ndvgrd,
     <                     rxydzr,rxydzzr)
              DO ip=1,nmzxy
                 dzmagmom(i,ip)= rxydzr(ip)
                ddzmagmom(i,ip)= rxydzzr(ip)
              ENDDO
            END DO
            DEALLOCATE ( rhtxyr,rxydzr,rxydzzr,dummy )
          END IF ! l_noco 

        ENDIF   ! igrd.GT.0

!       WRITE(6,'('' 9990nmzxy='',2i5)') nmzxy
        ALLOCATE ( rhdx(0:ifftd2-1,jspd),rhdy(0:ifftd2-1,jspd) )
        ALLOCATE ( rhdz(0:ifftd2-1,jspd),rhdxx(0:ifftd2-1,jspd) )
        ALLOCATE ( rhdyy(0:ifftd2-1,jspd),rhdzz(0:ifftd2-1,jspd) )
        ALLOCATE ( rhdyz(0:ifftd2-1,jspd),rhdzx(0:ifftd2-1,jspd) )
        ALLOCATE ( rhdxy(0:ifftd2-1,jspd),gggrd(0:ifftd2-1) )
        ALLOCATE ( agr(0:ifftd2-1),agru(0:ifftd2-1),agrd(0:ifftd2-1) )
        ALLOCATE ( g2r(0:ifftd2-1),g2ru(0:ifftd2-1),g2rd(0:ifftd2-1) )
        ALLOCATE ( gggr(0:ifftd2-1),gggru(0:ifftd2-1),gzgr(0:ifftd2-1) )
        ALLOCATE ( vxc(0:ifftd2-1,jspd),exc(0:ifftd2-1) )
        ALLOCATE ( cqpw(n2d-1),af2(0:ifftd2-1,jspd) )
        ALLOCATE ( fgxy(n2d-1),bf2(0:ifftd2-1) )

        rd = 0.0

        DO ip = 1,nmzxy
        ! loop over warping region

          IF (.not. l_noco) THEN
          ! Transform charge and magnetization to real-space.

            DO js=1,jspins
              CALL fft2d(
     >               k1d,k2d,n2d,
     <               af2(0,js),bf2,
     >               rht(ip,ivac,js),0.,rhtxy(ip,1,ivac,js),
     >               nmzxyd,ng2,kimax2,+1,
     >               igfft2,pgfft2,nstr2)
            END DO

          ELSE

            DO i=0,9*k1d*k2d-1
              af2(i,1)= af2noco(i,ip,1) 
              af2(i,2)= af2noco(i,ip,2)
            END DO

          END IF
 
          IF (igrd > 0) THEN 
          ! calculate derivatives with respect to x,y in g-space 
          ! and transform them to real-space.  

            DO js = 1,jspins

              DO iq=1,ng2-1
                cqpw(iq)=ci*rhtxy(ip,iq,ivac,js)
              ENDDO

              rhti = 0.0                    ! d(rho)/dx is obtained 
              CALL fft2d(                   ! by a FFT of i*gx*rhtxy
     >                 k1d,k2d,n2d,         ! (rht is set to zero,
     <                 rhdx(0,js),bf2,      !  and gx is included in
     >                 zro,rhti,cqpw,       !  pgft2x) 
     >                 1,ng2,kimax2,+1,     ! rhdx =
     >                 igfft2,pgft2x,nstr2) ! dn/dx =  FFT(0,i*gx*rhtxy)

              rhti = 0.0
              CALL fft2d(                   ! dn/dy =  FFT(0,i*gy*rhtxy)
     >                 k1d,k2d,n2d,
     <                 rhdy(0,js),bf2,
     >                 zro,rhti,cqpw,
     >                 1,ng2,kimax2,+1,
     >                 igfft2,pgft2y,nstr2)

              rhti = 0.0
              CALL fft2d(                   ! dn/dz = FFT(rhtdz,rxydz)
     >                 k1d,k2d,n2d,
     <                 rhdz(0,js),bf2,
     >                 rhtdz(ip,js),rhti,rxydz(ip,1,js),
     >                 nmzxyd,ng2,kimax2,+1,
     >                 igfft2,pgfft2,nstr2)


              DO iq=1,ng2-1
                cqpw(iq)=-rhtxy(ip,iq,ivac,js)
              ENDDO

              rhti = 0.0
              CALL fft2d(                ! d2n/dx2 = FFT(0,-gx^2*rhtxy)
     >                 k1d,k2d,n2d,
     <                 rhdxx(0,js),bf2,
     >                 zro,rhti,cqpw,
     >                 1,ng2,kimax2,+1,
     >                 igfft2,pgft2xx,nstr2)

              rhti = 0.0
              CALL fft2d(                 ! d2n/dy2 = FFT(0,-gy^2*rhtxy)
     >                 k1d,k2d,n2d,
     <                 rhdyy(0,js),bf2,
     >                 zro,rhti,cqpw,
     >                 1,ng2,kimax2,+1,
     >                 igfft2,pgft2yy,nstr2)

              rhti = 0.0
              CALL fft2d(                 ! d2n/dz2 = FFT(rhtdzz,rxydzz)
     >                 k1d,k2d,n2d,
     <                 rhdzz(0,js),bf2,
     >                 rhtdzz(ip,js),rhti,rxydzz(ip,1,js),
     >                 nmzxyd,ng2,kimax2,+1,
     >                 igfft2,pgfft2,nstr2)


              DO iq=1,ng2-1
                cqpw(iq)=ci*rxydz(ip,iq,js)
              ENDDO

              rhti = 0.0
              CALL fft2d(                  ! d2n/dyz = FFT(0,i*gy*rxydz)
     >                 k1d,k2d,n2d,
     <                 rhdyz(0,js),bf2,
     >                 zro,rhti,cqpw,
     >                 1,ng2,kimax2,+1,
     >                 igfft2,pgft2y,nstr2)

              rhti = 0.0
              CALL fft2d(                  ! d2n/dzx = FFT(0,i*gx*rxydz)
     >                 k1d,k2d,n2d,
     <                 rhdzx(0,js),bf2,
     >                 zro,rhti,cqpw,
     >                 1,ng2,kimax2,+1,
     >                 igfft2,pgft2x,nstr2)

              DO iq=1,ng2-1
                cqpw(iq)=-rhtxy(ip,iq,ivac,js)
              ENDDO

              rhti = 0.0
              CALL fft2d(               ! d2n/dxy = FFT(0,-gx*gy*rhtxy)
     >                 k1d,k2d,n2d,
     <                 rhdxy(0,js),bf2,
     >                 zro,rhti,cqpw,
     >                 1,ng2,kimax2,+1,
     >                 igfft2,pgft2xy,nstr2)

            END DO ! js=1,jspins


            IF (l_noco) THEN
c ! In non-collinear calculations the derivatives of |m| are calculated
c ! in real-space. The derivatives of the charge density, that are 
c ! already calculated in g-space, will be used. 

              CALL grdrsvac(
     >               magmom(0,ip),bmat1,3*k1d,3*k2d,ndvgrd,
     <               dxmagmom,dymagmom) 
              DO i=0,9*k1d*k2d-1 
                chdens= rhdx(i,1)/2.+rhdx(i,2)/2.
                rhdx(i,1)= chdens + dxmagmom(i) 
                rhdx(i,2)= chdens - dxmagmom(i) 
                chdens= rhdy(i,1)/2.+rhdy(i,2)/2.
                rhdy(i,1)= chdens + dymagmom(i) 
                rhdy(i,2)= chdens - dymagmom(i) 
                chdens= rhdz(i,1)/2.+rhdz(i,2)/2.
                rhdz(i,1)= chdens + dzmagmom(i,ip) 
                rhdz(i,2)= chdens - dzmagmom(i,ip) 
              END DO  

              CALL grdrsvac(
     >               dxmagmom,bmat1,3*k1d,3*k2d,ndvgrd, 
     <               ddxmagmom(0,1),ddymagmom(0,1))
              CALL grdrsvac(
     >               dymagmom,bmat1,3*k1d,3*k2d,ndvgrd, 
     <               ddxmagmom(0,2),ddymagmom(0,2))
              DO i=0,9*k1d*k2d-1
                chdens= rhdxx(i,1)/2.+rhdxx(i,2)/2. 
                rhdxx(i,1)= chdens + ddxmagmom(i,1) 
                rhdxx(i,2)= chdens - ddxmagmom(i,1) 
                chdens= rhdyy(i,1)/2.+rhdyy(i,2)/2. 
                rhdyy(i,1)= chdens + ddymagmom(i,2) 
                rhdyy(i,2)= chdens - ddymagmom(i,2) 
                chdens= rhdxy(i,1)/2.+rhdxy(i,2)/2. 
                rhdxy(i,1)= chdens + 
     &           ( ddxmagmom(i,2) + ddymagmom(i,1) )/2.
                rhdxy(i,2)= chdens - 
     &           ( ddxmagmom(i,2) + ddymagmom(i,1) )/2.
              END DO
              CALL grdrsvac(
     >               dzmagmom(0,ip),bmat1,3*k1d,3*k2d,ndvgrd, 
     <               ddxmagmom(0,1),ddymagmom(0,1))
              DO i=0,9*k1d*k2d-1
                chdens= rhdzx(i,1)/2.+rhdzx(i,2)/2. 
                rhdzx(i,1)= chdens + ddxmagmom(i,1) 
                rhdzx(i,2)= chdens - ddxmagmom(i,1) 
                chdens= rhdyz(i,1)/2.+rhdyz(i,2)/2. 
                rhdyz(i,1)= chdens + ddymagmom(i,1) 
                rhdyz(i,2)= chdens - ddymagmom(i,1) 
                chdens= rhdzz(i,1)/2.+rhdzz(i,2)/2. 
                rhdzz(i,1)= chdens + ddzmagmom(i,ip) 
                rhdzz(i,2)= chdens - ddzmagmom(i,ip) 
              END DO  

            END IF ! l_noco  

          END IF ! igrd > 0 
c
c set minimal value of af2 to 1.0e-13
c
          CALL setmn(
     >               ifftd2,kimax2,jspins,
     X               af2)


c   calculate the quantities such as abs(grad(rho)),.. used in
cc  evaluating the gradient contributions to potential and energy.

         rd = z1 + delz*(ip-1)

         if(odi%d1)then
            CALL od_mkgxyz3(
     >           igrd,ifftd2,jspd,ifftd2,jspins,
     >           af2,rd,rhdx,rhdy,rhdz,rhdxx,rhdyy,
     >           rhdzz,rhdyz,rhdzx,rhdxy,
     <           agr,agru,agrd,g2r,g2ru,g2rd,
     <           gggr,gggru,gggrd,gzgr)
         else
            CALL mkgxyz3(
     >           igrd,ifftd2,jspd,ifftd2,jspins,af2,rhdx,rhdy,
     >           rhdz,rhdxx,rhdyy,rhdzz,rhdyz,rhdzx,rhdxy,
     <           agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,
     <           gzgr)
         endif

c     rhmnv: rho_minimum_vacuum.

          rhmnv=10.e+10
          DO js=1,jspins
            DO i=0,kimax2
              af2(i,js)=max(af2(i,js),d_15)
              rhmnv=min(rhmnv,af2(i,js))
            ENDDO
          ENDDO

          IF (rhmnv.lt.rhmn) THEN
            rhmn=rhmnv
            ichsmrg=3
          ENDIF

          IF (rhmn.LT.chng) THEN
            WRITE(6,'(/'' rhmn.lt.chng. rhmn,chng='',2d9.2)') rhmn,chng
c            STOP 'vvacxcg: rhmn.lt.chng'
          ENDIF

c         calculate the exchange-correlation potential in  real space
c

          CALL vxcallg(
     >                 icorr,lwbc,jspins,ifftd2,nt,af2,agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 vxc,
     >                 idsprs,isprsv,sprsv)


          CALL dset(ifftd2,0.0,bf2,1)

          DO js = 1,jspins
c
c           ----> 2-d back fft to g space
c
            CALL dset(ifftd2,0.0,bf2,1)
            CALL fft2d(
     >                 k1d,k2d,n2d,
     >                 vxc(0,js),bf2,
     <                 fgz,rhti,fgxy,
     >                 1,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)

c            ----> and add vxc to coulomb potential
c                  the g||.eq.zero component is added to vz
c
             vz(ip,ivac,js) = fgz + vz(ip,ivac,js)
c
c            the g||.ne.zero components are added to vxy
c
             DO irec2 = 1,ng2-1
               vxy(ip,irec2,ivac,js)=vxy(ip,irec2,ivac,js)+fgxy(irec2)
             ENDDO

          END DO   


c         calculate the exchange-correlation energy density in  real space
c
          IF (total) THEN

            CALL excallg(
     >                   icorr,lwbc,jspins,ifftd2,nt,af2,agr,agru,agrd,
     <                   g2r,g2ru,g2rd, gggr,gggru,gggrd,gzgr,
     <                   exc,
     >                   idsprs,isprsv,sprsv)

c           ----> 2-d back fft to g space
c
            CALL dset(ifftd2,0.0,bf2,1)
            CALL fft2d(
     >                 k1d,k2d,n2d,
     >                 exc,bf2,
     <                 excz(ip,ivac),rhti,excxy(ip,1,ivac),
     >                 nmzxyd,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)

          ENDIF

        END DO ! ip=1,nmzxy 
        DEALLOCATE ( rhdx,rhdy,rhdz,rhdxx,rhdyy,rhdzz )
        DEALLOCATE ( cqpw,fgxy,     rhdyz,rhdzx,rhdxy )

        IF (l_noco) THEN
          DEALLOCATE ( dzmagmom,ddzmagmom,dxmagmom,af2noco )
          DEALLOCATE ( dymagmom,ddxmagmom,ddymagmom )
        END IF  

        ! now treat the non-warping region 

        nmzdiff = nmz - nmzxy
c       WRITE(6,'(/'' 9992excz''/(8f15.7))') (excz(ip,1),ip=1,nmz)
        WRITE(6,'(/'' 9992nmzdiff='',i5)') nmzdiff

        ! The non-warping region runs from nmzxy+1 to nmz.
        ! The values from nmz0 to nmzxy are taken into account in order
        ! to get the real-space derivative smooth around nmzxy+1. 
        nmz0= nmzxy+1+(ndvgrd/2)-ndvgrd
        IF (nmz0 <= 0) THEN ! usually nmzxy>ndvgrd 
          nmz0= 1
        END IF 

        ALLOCATE ( rhtz(nmzd,jspd) )

        DO ip=nmz0,nmz 
          IF (.not. l_noco) THEN
            DO js=1,jspins 
              rhtz(ip,js)= rht(ip,ivac,js)
            END DO
          ELSE
            af2(0,1) = rht(ip,ivac,1)
            af2(0,2) = rht(ip,ivac,2)
            mx(0)= REAL(cdomvz(ip,ivac))
            my(0)= AIMAG(cdomvz(ip,ivac))
            chdens= (af2(0,1)+af2(0,2))/2.
            magmom(0,1)= mx(0)**2 + my(0)**2 +
     &                   ((af2(0,1)-af2(0,2))/2.)**2
            magmom(0,1)= SQRT(magmom(0,1))
            rhtz(ip,1)= chdens + magmom(0,1)
            rhtz(ip,2)= chdens - magmom(0,1) 
          END IF 
        END DO 

        IF (l_noco) THEN 
          DEALLOCATE ( magmom,mx,my ) 
          ALLOCATE ( dummy(nmz) )
          DO js=1,jspins
            CALL grdchlh(
     >                   0,1,nmz-nmz0+1,delz,dummy,rhtz(nmz0,js),ndvgrd,
     <                   rhtdz(nmz0,js),rhtdzz(nmz0,js))
          END DO
          DEALLOCATE ( dummy )

        END IF 

c       calculate the quantities such as abs(grad(rho)),.. used in
cc      evaluating the gradient contributions to potential and
cc      energy.

        agr(:)=0.0 ; agru(:)=0.0 ; agrd(:)=0.0 ; g2r(:)=0.0
        g2ru(:)=0.0 ; g2rd(:)=0.0 ; gggr(:)=0.0 ; gggru(:)=0.0
        gggrd(:)=0.0 ; gzgr(:)=0.0

        IF (igrd.gt.0)  THEN
          if(odi%d1)then
             CALL od_mkgz(
     >            z1,nmzxy,delz,
     >            nmzdiff,jspins,
     >            rhtz(nmzxy+1,1),rhtz(nmzxy+1,jspins),
     >            rhtdz(nmzxy+1,1), rhtdz(nmzxy+1,jspins),
     >            rhtdzz(nmzxy+1,1),rhtdzz(nmzxy+1,jspins),
     <            agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,
     <            gzgr)

          else
             CALL mkgz(nmzdiff,jspins,
     >            rhtz(nmzxy+1,1),rhtz(nmzxy+1,jspins),
     >            rhtdz(nmzxy+1,1),rhtdz(nmzxy+1,jspins),
     >            rhtdzz(nmzxy+1,1),rhtdzz(nmzxy+1,jspins),
     <            agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,
     <            gzgr)
          endif
        ENDIF

c       calculate vxc for z now beyond warping region

        ALLOCATE ( rhtxc(nmzd,jspd),vxcz(nmzd,jspd) )

        DO js=1,jspins
          DO i=0,ifftd2-1
            af2(i,js)=0.0
          ENDDO
          DO ip=nmzxy+1,nmz
            rhtxc(ip-nmzxy,js) = max(rhtz(ip,js),d_15) !+gb
            rhmnv=min(rhmnv,rhtxc(ip-nmzxy,js))
          ENDDO
        ENDDO

        IF(rhmnv.lt.rhmn) THEN
          rhmn=rhmnv
          ichsmrg=4
        ENDIF

        IF (rhmn.lt.chng) THEN
          WRITE (6,'('' rhmn.lt.chng. rhmn,chng='',2d9.2)') rhmn,chng
c          STOP 'vvacxcg: rhmn.lt.chng'
        ENDIF

        CALL vxcallg(
     >               icorr,lwbc,jspins,nmzd,nmzdiff,rhtxc,agr,agru,
     >               agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <               vxcz,
     >               idsprs,isprsv,sprsv)

        DO js = 1,jspins
           DO ip = nmzxy + 1,nmz
              vz(ip,ivac,js) = vz(ip,ivac,js) + vxcz(ip-nmzxy,js)
           ENDDO
        ENDDO

c
        WRITE (6,fmt=8020) ivac, (vz(nmz,ivac,js),js=1,jspins)
        WRITE(16,fmt=8020) ivac, (vz(nmz,ivac,js),js=1,jspins)
 8020   FORMAT(/,5x,'vacuum zero for vacuum',i3,' = ',2f14.10)
c
c     calculate the ex-corr. energy density now beyond warping region
c
         IF (total) THEN
            CALL excallg(
     >                   icorr,lwbc,jspins,nmzd,nmzdiff,rhtxc,agr,agru,
     >                   agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                   excz(nmzxy+1,ivac),
     >                   idsprs,isprsv,sprsv)
         ENDIF

        DEALLOCATE ( af2,bf2,agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru )
        DEALLOCATE ( gggrd,gzgr,vxc,exc,vxcz,rhtz,rhtxc )

      ENDDO    ! loop over vacua (ivac)
      DEALLOCATE ( rhtdz,rhtdzz,rxydz,rxydzz )

      idsprs=idsprssv

      END SUBROUTINE vvacxcg
      END MODULE m_vvacxcg
