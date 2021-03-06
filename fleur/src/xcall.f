      MODULE m_xcall
!***********************************************************************
!
!  vxcall ... driver subroutine for the XC potential and
!  excall ... (only in case of total = .true.)  the XC energy density
!
!  INPUT :  electron density rh()                        (hartree units)
!  OUTPUT:  effective exch-correlation potential vxc()   (hartree units)
!                               energy density   exc()
!
!   icorr selects various exchange correlation potentials:
!   icorr =            0 means X-alpha method
!                      1 means Wigner correlation
!                      2 means Moruzzi,Janak,Williams correlat
!                      3 means von Barth,Hedin correlation
!                      4 means Vosko,Wilk,Nusair correlation
!                      5 means Perdew,Zunger correlation
!
!   krla=1 : relativistic correction according to
!            A. H. MacDonald and S. H. Vosko, J. Phys. C12, 2977 (1979)
!
!   based on a subroutine from S.Bluegel,  R.Pentcheva, 22.01.96
!
!***********************************************************************
      IMPLICIT NONE
      INTEGER, PRIVATE :: ir
      REAL, PARAMETER, PRIVATE :: hrtr_half = 0.5

      CONTAINS
!***********************************************************************
      SUBROUTINE vxcall(
     >                  iofile,icorr,krla,jspins,
     >                  mgrid,ngrid,rh,
     <                  vxc)
!***********************************************************************

      USE m_xcxal, ONLY : vxcxal
      USE m_xcwgn, ONLY : vxcwgn
      USE m_xcbh,  ONLY : vxcbh
      USE m_xcvwn, ONLY : vxcvwn
      USE m_xcpz,  ONLY : vxcpz
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: iofile              ! file number for read and write
      INTEGER, INTENT (IN) :: icorr,krla,jspins   ! run mode parameters
      INTEGER, INTENT (IN) :: ngrid,mgrid         ! mesh,number of mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: vxc(mgrid,jspins)     ! xc energy  potential
!
!--> Determine exchange correlation energy density
!
      IF (icorr.EQ.0)  THEN   ! X-alpha method
         CALL vxcxal(
     >               krla,jspins,
     >               mgrid,ngrid,rh,
     <               vxc)

      ELSEIF (icorr.EQ.1) THEN    ! Wigner interpolation formula
         CALL vxcwgn(
     >               krla,jspins,
     >               mgrid,ngrid,rh,
     <               vxc)
      ELSEIF ((icorr.EQ.2).or.(icorr.EQ.3)) THEN ! von Barth,Hedin correlation
        CALL vxcbh(
     >             iofile,icorr,krla,jspins,
     >             mgrid,ngrid,rh,
     <             vxc)

      ELSEIF (icorr.EQ.4) THEN     ! Vosko,Wilk,Nusair correlation
        CALL vxcvwn(
     >              iofile,krla,jspins,
     >              mgrid,ngrid,rh,
     <              vxc)
      ELSEIF (icorr.EQ.5) THEN     ! Perdew,Zunger correlation
        CALL vxcpz(                     
     >             iofile,krla,jspins,
     >             mgrid,ngrid,rh,
     <             vxc)
      ELSE
         WRITE (iofile,FMT=9000) icorr
         STOP 'vxcall'
 9000    FORMAT (13x,'set key for exchange-correlation potential',i2)
      ENDIF
!
!---> convert to hartree units
!
      IF (jspins.EQ.2) THEN
         DO ir = 1,ngrid
            vxc(ir,1)      = hrtr_half * vxc(ir,1)
            vxc(ir,jspins) = hrtr_half * vxc(ir,jspins)
         ENDDO
      ELSE IF (jspins.EQ.1) THEN
         DO ir = 1,ngrid
            vxc(ir,1) = hrtr_half * vxc(ir,1)
         ENDDO
      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
         STOP 'vxcall'
      END IF

      RETURN
      END SUBROUTINE vxcall

!***********************************************************************
      SUBROUTINE  excall(
     >                   iofile,icorr,krla,jspins,
     >                   mgrid,ngrid,rh,
     <                   exc)
!***********************************************************************

      USE m_xcxal, ONLY : excxal
      USE m_xcwgn, ONLY : excwgn
      USE m_xcbh,  ONLY : excbh
      USE m_xcvwn, ONLY : excvwn
      USE m_xcpz,  ONLY : excpz
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: iofile              ! file number for read and write
      INTEGER, INTENT (IN) :: icorr,krla,jspins   ! run mode parameters
      INTEGER, INTENT (IN) :: ngrid,mgrid         ! mesh,number of mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy  density
!
!--> Determine exchange correlation energy density
!
      IF (icorr.EQ.0)  THEN   ! X-alpha method
         CALL excxal(
     >               iofile,krla,jspins,
     >               mgrid,ngrid,rh,
     <               exc)
 
      ELSEIF (icorr.EQ.1) THEN    ! Wigner interpolation formula
         CALL excwgn(
     >               iofile,krla,jspins,
     >               mgrid,ngrid,rh,
     <               exc)
      ELSEIF ((icorr.EQ.2).or.(icorr.EQ.3)) THEN ! von Barth,Hedin correlation
        CALL excbh(
     >             iofile,icorr,krla,jspins,
     >             mgrid,ngrid,rh,
     <             exc)

      ELSEIF (icorr.EQ.4) THEN     ! Vosko,Wilk,Nusair correlation
        CALL excvwn(
     >              iofile,krla,jspins,
     >              mgrid,ngrid,rh,
     <              exc)
      ELSEIF (icorr.EQ.5) THEN     ! Perdew,Zunger correlation
         CALL excpz(                      
     >              iofile,krla,jspins,
     >              mgrid,ngrid,rh,
     <              exc)
      ELSE
         WRITE (iofile,FMT=9001)
         STOP 'excall'
 9001    FORMAT (13x,'set key for exchange-correlation potential')
      ENDIF
!
!---> convert to hartree units
!
      DO ir = 1,ngrid
         exc(ir) = hrtr_half * exc(ir)
      ENDDO

      RETURN
      END SUBROUTINE excall
      END MODULE m_xcall
