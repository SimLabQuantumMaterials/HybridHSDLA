      MODULE m_visxc
c     ******************************************************
c     subroutine generates the exchange-correlation potential
c     in the interstitial region    c.l.fu
c     ******************************************************
      CONTAINS
      SUBROUTINE visxc(
     >                 ifftd,k1d,k2d,k3d,n3d,jspd,
     >                 qpw,cdom,l_noco,
     >                 kimax,igfft,pgfft,ufft,
     >                 icorr,total,krla,jspins,ng3,nstr,
     X                 vpw,vpw_w,
     <                 excpw)

c     ******************************************************
c     instead of visxcor.f: the different exchange-correlation 
c     potentials defined through the key icorr are called through 
c     the driver subroutine vxcall.f,for the energy density - excall
c     subroutines vectorized
c     in case of TOTAL = .TRUE. calculates the ex.-corr. energy density
c     ** r.pentcheva 08.05.96
c     ********************************************************************

      USE m_xcall, ONLY : vxcall,excall
      USE m_fft3d
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      LOGICAL, INTENT (IN) :: l_noco 
      INTEGER, INTENT (IN) :: ifftd,k1d,k2d,k3d,n3d,jspd
      INTEGER, INTENT (IN) :: kimax,icorr,krla,jspins,ng3
      LOGICAL, INTENT (IN) :: total
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      INTEGER, INTENT (IN) :: nstr(n3d)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
      COMPLEX, INTENT (IN) :: qpw(n3d,jspd),cdom(n3d)  
      COMPLEX, INTENT (OUT) :: excpw(n3d)
      COMPLEX, INTENT (INOUT) ::vpw(n3d,jspd),vpw_w(n3d,jspd)
C     ..
C     .. Local Scalars ..
      INTEGER i,k,js,nt,ivec(n3d)
      REAL    chdens,magmom
C     ..
C     .. Local Arrays ..
      COMPLEX fg3(n3d)
      REAL, ALLOCATABLE :: mx(:),my(:)
      REAL, ALLOCATABLE :: exc(:),vcon(:),vxc(:,:)
      REAL, ALLOCATABLE :: af3(:,:),bf3(:)
c
c     ---> allocate arrays
c
      ALLOCATE ( exc(0:ifftd-1),vcon(0:ifftd-1),vxc(0:ifftd-1,jspd),
     +    af3(0:ifftd-1,jspd),bf3(0:ifftd-1) )

c
c     transform charge density to real space
c
      DO js = 1,jspins
         CALL fft3d(
     <              af3(0,js),bf3,
     >              qpw(1,js),
     >              k1d,k2d,k3d,ng3,kimax,+1,
     >              igfft,pgfft,nstr)
      ENDDO

      IF (l_noco) THEN 

        ALLOCATE (mx(0:ifftd-1),my(0:ifftd-1))

        CALL fft3d(
     <             mx,my,
     >             cdom,
     >             k1d,k2d,k3d,ng3,kimax,+1,
     >             igfft,pgfft,nstr)
        DO i=0,27*k1d*k2d*k3d-1
          chdens= (af3(i,1)+af3(i,2))/2.
          magmom= mx(i)**2 + my(i)**2 + ((af3(i,1)-af3(i,2))/2.)**2
          magmom= SQRT(magmom)
          af3(i,1)= chdens + magmom
          af3(i,2)= chdens - magmom
        END DO

        DEALLOCATE (mx,my)

      END IF

c
c     calculate the exchange-correlation potential in  real space
c
       nt=ifftd
       CALL vxcall
     >            (6,icorr,krla,jspins,
     >             ifftd,nt,af3,
     <             vxc)

c
c     ---> back fft to g space and add to coulomb potential for file nrp
c

      IF (total) THEN

         ivec(1:ng3)=1
         DO js = 1,jspins

           DO i=0,ifftd-1
             vcon(i)=ufft(i)*vxc(i,js)
             bf3(i)=0.0
           ENDDO
           CALL fft3d(
     >                vcon,bf3,
     <                fg3,
     >                k1d,k2d,k3d,ng3,kimax,-1,
     >                igfft,pgfft,ivec)

           DO k = 1,ng3
              vpw_w(k,js) = vpw_w(k,js) + fg3(k)
           ENDDO   

         ENDDO
c
c     calculate the ex.-cor energy density in real space
c
          CALL excall
     >               (6,icorr,krla,jspins,
     >                ifftd,nt,af3,
     <                exc)

         DO i=0,ifftd-1
           vcon(i)=ufft(i)*exc(i)
           bf3(i)=0.0
         ENDDO
c
c         ---> back fft to g space
c
         CALL fft3d(
     >              vcon,bf3,
     <              excpw,
     >              k1d,k2d,k3d,ng3,kimax,-1,
     >              igfft,pgfft,ivec)
c
      END IF ! total

      DO js = 1,jspins
         bf3(0:ifftd-1)=0.0
         CALL fft3d(
     >              vxc(0,js),bf3,
     <              fg3,
     >              k1d,k2d,k3d,ng3,kimax,-1,
     >              igfft,pgfft,nstr)
         DO k = 1,ng3
            vpw(k,js) = vpw(k,js) + fg3(k)
         ENDDO   
      ENDDO

      DEALLOCATE ( exc,vcon,vxc,af3,bf3 )

      END SUBROUTINE visxc
      END MODULE m_visxc
