      MODULE m_xcbh 
!-----------------------------------------------------------------------
!     Called in case of icorr=2,3 : spin-polarized exchange-correlation 
!       potential of U. von Barth and L. Hedin, J.Phys.C5,1629 (1972)
!       icorr = 2: parametrization of Moruzzi,Janak,Williams
!       icorr = 3: parametrization of von Barth and Hedin
!
!     krla=1: Relativistic correction of exchange energy and potential 
!             related to Dirac kinetic energy, according to:
!             A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)
!
!     be careful: calculation in rydberg!
!
!     vxcbh calculates the XC-potential and
!     excbh calculates the XC-energy
!
!     based on a subroutine by S. Bluegel;   r.pentcheva 22.01.96
!-----------------------------------------------------------------------

      USE m_constants, ONLY : pimach
      USE m_relcor
      IMPLICIT NONE

      REAL, PARAMETER, PRIVATE :: ff  = 3.847322101863  ! 1 / ( 2^(1/3) - 1 )
      REAL, PARAMETER, PRIVATE :: cvx = 1.221774115422  ! 2 * ( 3/(2*pi) )^(2/3)
      REAL, PARAMETER, PRIVATE :: cpmjw = 0.045  , cfmjw = 0.0225
      REAL, PARAMETER, PRIVATE :: cpvbh = 0.0504 , cfvbh = 0.0254
      REAL, PARAMETER, PRIVATE :: rpmjw = 21.0 , rfmjw = 52.916684096
      REAL, PARAMETER, PRIVATE :: rpvbh = 30.0 , rfvbh = 75.0
      REAL, PARAMETER, PRIVATE :: d_15 = 1.e-15
      REAL, PARAMETER, PRIVATE :: one = 1.0 , three = 3.0 , four = 4.0
      REAL, PARAMETER, PRIVATE :: half = 0.5 , thrd = one/three
      REAL, PARAMETER, PRIVATE :: hfthrd = 0.79370052705 ! 2^(-1/3)
      REAL, PARAMETER, PRIVATE :: thrhalf = three * half
      REAL, PARAMETER, PRIVATE :: fothrd = four * thrd , two = 2.0

      REAL, PRIVATE :: rho, rh1, rh2 ! total, spin up & spin down charge density
      REAL, PRIVATE :: pi, x, y, cp, cf, rp, rf, rs, ecprs, ecfrs
      INTEGER, PRIVATE :: ir

      CONTAINS
!************************************************************************
      SUBROUTINE vxcbh
     >                (iofile,icorr,krla,jspins,
     >                 mgrid,ngrid,rh,
     <                 vxc)
!************************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: icorr,krla  !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: vxc(mgrid,jspins)     ! xc potential
!
!     .. Local Scalars ..
      REAL txthrd,tythrd,muxp,mucp,mucf,ecfmp,tauc,mucnm
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: psi(:)       ! relativistic exchange potential corr.
!
!---- Intrinsic Functions
      INTRINSIC alog,max

      pi = pimach() 
!
!-----> evaluate relativistic corrections for exchange
! 
      ALLOCATE ( psi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.true.,rh,
     <            psi)
!
!-----> select exchange correlation potential
!
      IF (icorr.EQ.2) THEN
         cp = cpmjw ; cf = cfmjw
         rp = rpmjw ; rf = rfmjw
      ELSEIF (icorr.EQ.3) THEN
         cp = cpvbh ; cf = cfvbh
         rp = rpvbh ; rf = rfvbh
      ELSE
         WRITE (iofile,FMT=2000)
         STOP 'vxcbh'
      END IF
 2000 FORMAT (13x,'set key for exchange-correlation potential')
!
!-----> calculate exchange correlation potential
!
      IF ( jspins .EQ. 2) THEN               ! spinpolarized calculation

        DO ir = 1,ngrid                        ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1))
          rh2 = max(d_15,rh(ir,jspins))
          rho = rh1 + rh2
          x = rh1/rho
          y = rh2/rho
          txthrd = (2*x)**thrd
          tythrd = (2*y)**thrd
          rs= (four*pi*rho/three)**thrd
          rs = 1/rs

          ecprs = -cp*fc(rs/rp)            ! calculate correlation energy
          ecfrs = -cf*fc(rs/rf)            ! p : paramagnetic, f : ferromagnetic
                                           ! x : exchange,     c : correlation
          muxp = -psi(ir)* (cvx/rs)        ! paramagnetic exchange potential 
                                           !       (psi contains rel. corr.)
          mucp = -cp*alog(one+rp/rs)       ! calculate correlation potential
          mucf = -cf*alog(one+rf/rs)
          ecfmp = fothrd * (ecfrs-ecprs)
          tauc = mucf - mucp - ecfmp
          mucnm = mucp + tauc*fex(x) - ff*ecfmp

          vxc(ir,1)      = mucnm + (muxp+ff*ecfmp)*txthrd  ! collect correlation 
          vxc(ir,jspins) = mucnm + (muxp+ff*ecfmp)*tythrd  ! and exchange parts
        ENDDO

      ELSEIF ( jspins .EQ. 1 ) THEN        ! non - spinpolarized calculation

        DO ir = 1, ngrid                   ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1))
          rs = (four*pi*rh1/three)**thrd
          rs = 1/rs 
          muxp = -psi(ir) * (cvx/rs)       ! paramagnetic exchange potential 
                                           !       (psi contains rel. corr.)
          mucp = -cp* alog(one+rp/rs)      ! calculate correlation potential
          vxc(ir,1)     = mucp + muxp      ! collect correlation & exchange part
        ENDDO

      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
         STOP 'vxcbh'
      ENDIF

      DEALLOCATE (psi)
      RETURN   

      END SUBROUTINE vxcbh
C***********************************************************************
      SUBROUTINE excbh
     >                (iofile,icorr,krla,jspins,
     >                 mgrid,ngrid,rh,
     <                 exc)
C***********************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: icorr,krla  !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy
!
!     .. Local Scalars ..
      REAL thfpi,thrquart,exprs,exfrs,excprs,excfrs
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: phi(:)       ! relativistic exchange energy correct.
!
!-----> Intrinsic Functions
      INTRINSIC alog,max
!
      pi = pimach()
      thrquart = 0.75
      thfpi = thrquart/pi

      ALLOCATE ( phi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.false.,rh,
     <            phi)
!
!-----> select exchange correlation potential
!
      IF (icorr.EQ.2) THEN
         cp = cpmjw ; cf = cfmjw
         rp = rpmjw ; rf = rfmjw
      ELSEIF (icorr.EQ.3) THEN
         cp = cpvbh ; cf = cfvbh
         rp = rpvbh ; rf = rfvbh
      ELSE
         WRITE (iofile,FMT=2000)
         STOP 'excbh'
      END IF
 2000 FORMAT (13x,'set key for exchange-correlation potential')

      IF ( jspins .EQ. 2) THEN       ! spinpolarized calculation

        DO  ir = 1,ngrid                  ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1))
          rh2 = max(d_15,rh(ir,jspins))
          rho = rh1 + rh2
          x = rh1/rho
          rs= (thfpi/rho)**thrd

          exprs = -phi(ir)*thrquart*cvx/rs    ! first exchange energy 
          exfrs = (2.0**thrd)*exprs           ! phi contains rel. corr.

          ecprs = -cp*fc(rs/rp)               ! calculate correlation energy
          ecfrs = -cf*fc(rs/rf)               ! p: paramagnetic, f: ferromagn.

          excprs = exprs + ecprs              ! now add correlation energy
          excfrs = exfrs + ecfrs

          exc(ir) = excprs + (excfrs-excprs)*fex(x) ! collect all terms
        ENDDO 

      ELSEIF ( jspins .EQ. 1 ) THEN  ! non - spinpolarized calculation

        DO ir = 1,ngrid                    ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1)) 
          rs = (thfpi/rh1)**thrd
          exprs = -phi(ir)*thrquart*cvx/rs ! exchange energy ; phi contains 
                                           ! relativistic correctionS
          ecprs = -cp*fc(rs/rp)            ! calculate correlation energy
          exc(ir) = exprs + ecprs          ! add correlation energy
        ENDDO

       ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
         STOP 'excbh'
       ENDIF

      DEALLOCATE (phi)
      RETURN

      END SUBROUTINE excbh
!--------------------------------------------------------------------
      REAL FUNCTION fc(x)
        REAL x
        fc = (one+(x)*(x)*(x))*alog(one+one/(x))
     +        + half*(x) - (x)*(x) - thrd
      END  FUNCTION fc
      REAL FUNCTION fex(x)
        REAL x
        fex = ff/hfthrd*((x)**fothrd +(1-(x))**fothrd - hfthrd)
      END FUNCTION fex
!--------------------------------------------------------------------

      END MODULE m_xcbh
