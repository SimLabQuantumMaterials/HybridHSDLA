      MODULE m_xcwgn
!-----------------------------------------------------------------------
!     Called in case of icorr=1 : non-spinpolarized exchange-correlation 
!                                     from wigners interpolation formula.
!
!     krla=1: Relativistic correction of exchange energy and potential 
!             related to Dirac kinetic energy, according to:
!             A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)
!
!     be careful: calculation in rydberg!
!
!     vxc   calculates the XC-potential and
!     exc   calculates the XC-energy
!
!     based on a subroutine by S. Bluegel;   r.pentcheva 22.01.96
!-----------------------------------------------------------------------

      USE m_constants, ONLY : pimach
      USE m_relcor
      IMPLICIT NONE


      REAL, PARAMETER, PRIVATE :: one = 1.0 , three = 3.0 , four = 4.0
      REAL, PARAMETER, PRIVATE :: thrd = one/three , d_15 = 1.e-15
      REAL, PARAMETER, PRIVATE :: cex = 0.91633058742  ! 3/2 * ( 3/(2*pi) )^(2/3)

      REAL, PRIVATE :: pi, rs, rho, thfpi, exp, ecp
      INTEGER, PRIVATE :: ir

      CONTAINS
!************************************************************************
      SUBROUTINE vxcwgn(
     >                  krla,jspins,
     >                  mgrid,ngrid,rh,
     <                  vxc)
!************************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  run mode parameters
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: vxc(mgrid,jspins)     ! xc potential
!
!     .. Local Scalars ..
      REAL fothrd,vxp,vcp
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: psi(:)       ! relativistic exchange potential corr.
!
!-----s Intrinsic Functions
      INTRINSIC alog,max

      pi = pimach() 
      thfpi  = three / ( four * pi )
      fothrd = four/three
!
!-----> evaluate relativistic corrections for exchange
!
      ALLOCATE ( psi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.true.,rh,
     <            psi)

      DO ir = 1,ngrid
         IF (jspins.EQ.1) THEN
           rho = max(d_15,rh(ir,1))
         ELSE
           rho = max(d_15,rh(ir,1)) + max(d_15,rh(ir,jspins))
         ENDIF
         rs  = (thfpi/rho)**thrd
         exp = - cex / rs              ! exchange energy = -0.9163306/rs
         ecp = - 0.88 / (rs + 7.8)     ! correlation energy = -0.88/(rs + 7.8)

         vxp = fothrd * exp            ! exchange potential = 4/3 exp
         vcp = fothrd * ecp + 2.288 / ( rs + 7.8 ) ** 2

         vxc(ir,1) = vcp + vxp * psi(ir)
      ENDDO
      IF ( jspins .EQ. 2 ) THEN        ! spinpolarized calculation
         DO ir = 1,ngrid
            vxc(ir,jspins) = vxc(ir,1)
         ENDDO
      ENDIF

      DEALLOCATE (psi)
      RETURN

      END SUBROUTINE vxcwgn
C***********************************************************************
      SUBROUTINE excwgn(
     >                  iofile,krla,jspins,
     >                  mgrid,ngrid,rh,
     <                  exc)
C***********************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: phi(:)       ! relativistic exchange energy correct.
!
!-----> Intrinsic Functions
      INTRINSIC alog,max
!
      pi = pimach()
      thfpi  = three / ( four * pi )

      ALLOCATE ( phi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.false.,rh,
     <            phi)

      DO ir = 1,ngrid
         IF (jspins.EQ.1) THEN
           rho = max(d_15,rh(ir,1))
         ELSE
           rho = max(d_15,rh(ir,1)) + max(d_15,rh(ir,jspins))
         ENDIF
         rs  = (thfpi/rho)**thrd
         exp = - cex / rs              ! exchange energy = -0.9163306/rs
         ecp = - 0.88 / (rs + 7.8)     ! correlation energy = -0.88/(rs + 7.8)

         exc(ir) = ecp + exp * phi(ir)
      ENDDO
      IF (jspins.EQ.2) THEN
         WRITE (iofile,'('' WARNING: Wigner correlation !'',
     +   ''only applicable for non-spinpolarized calculations'')')
      ENDIF

      DEALLOCATE (phi)
      RETURN

      END SUBROUTINE excwgn
      END MODULE m_xcwgn
