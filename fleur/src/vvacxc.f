      MODULE m_vvacxc
c     ********************************************************************
c     calculates 2-dim star function coefficients of exchange-correlation*
c     potential in the vacuum regions  and adds them to the corresponding*
c     coeffs of the coulomb potential            c.l.fu, r.podloucky     *
c     ********************************************************************
      CONTAINS
      SUBROUTINE vvacxc(
     >                  k1d,k2d,nmzd,nmzxyd,n2d,jspd,ifftd2,
     >                  icorr,total,krla,nmzxy,jspins,ng2,nmz,nstr2,
     >                  rhtxy,rht,cdomvxy,cdomvz,l_noco,
     >                  kimax2,igfft2,pgfft2,nvac,
     X                  vxy,vz,
     <                  excxy,excz)

c     ********************************************************************
c     instead of vvacxcor.f: the different exchange-correlation 
c     potentials defined through the key icorr are called through 
c     the driver subroutine vxcall.f, subroutines vectorized
c     in case of TOTAL = .TRUE. calculates the ex.-corr. energy
c     density through the driver subroutine excall.f
c     ** r.pentcheva 08.05.96
c     ********************************************************************

      USE m_xcall, ONLY : vxcall,excall
      USE m_fft2d
      USE m_set, ONLY : dset

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,nmzd,nmzxyd,n2d,jspd,nvac,ifftd2
      INTEGER, INTENT (IN) :: icorr,krla,kimax2,nmzxy,jspins,ng2,nmz
      LOGICAL, INTENT (IN) :: l_noco,total
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nstr2(n2d)
      INTEGER, INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL,    INTENT (INOUT) :: rht(nmzd,2,jspd)
      COMPLEX, INTENT (INOUT) :: rhtxy(nmzxyd,n2d-1,2,jspd)
      COMPLEX, INTENT (INOUT) :: cdomvz(nmzd,2) 
      COMPLEX, INTENT (INOUT) :: cdomvxy(nmzxyd,n2d-1,2) 
      REAL,    INTENT (OUT) :: excz(nmzd,2)
      COMPLEX, INTENT (OUT) :: excxy(nmzxyd,n2d-1,2)
      REAL,    INTENT (INOUT) :: vz(nmzd,2,jspd)
      COMPLEX, INTENT (INOUT) :: vxy(nmzxyd,n2d-1,2,jspd) 
C     ..
C     .. Local Scalars ..
      INTEGER :: k,js,nt,irec2,nmzdiff,ivac,ip,i 
      REAL    :: rhti
      REAL    :: chdens,magmom 
C     ..
C     .. Local Arrays ..
      COMPLEX :: fgxy(n2d-1)
      REAL    :: af2(0:ifftd2-1,jspd),bf2(0:ifftd2-1),fgz
      REAL,ALLOCATABLE :: mx(:),my(:) 
c     warping region
      REAL    :: vxc(0:ifftd2-1,jspd),exc(0:ifftd2-1)
c     beyond warping region
      REAL    :: vxcz(nmzd,jspd)

      IF (l_noco) THEN
        ALLOCATE (mx(0:ifftd2-1),my(0:ifftd2-1)) 
      END IF 

      nt = ifftd2
      rhti = 0.
c
c     the charge density in vacuum is expanded in 2-dim stars on a mesh 
c     in z-direction . The G||.ne.zero-components expand from 1 to nmzxy
c     the G||.eq.zero-components expand from 1 to nmz
c     first we calculate vxc in the warping region
c
      DO 150 ivac = 1,nvac
        DO 110 ip = 1,nmzxy
c
c         transform charge density to real space: 2-dim FFT
c
          DO js = 1,jspins
            CALL fft2d(
     >               k1d,k2d,n2d,
     <               af2(0,js),bf2,
     >               rht(ip,ivac,js),rhti,rhtxy(ip,1,ivac,js),
     >               nmzxyd,ng2,kimax2,+1,
     >               igfft2,pgfft2,nstr2)
          END DO

          IF (l_noco) THEN 

            CALL fft2d(
     >               k1d,k2d,n2d,
     <               mx,my, 
     >               REAL(cdomvz(ip,ivac)),AIMAG(cdomvz(ip,ivac)),
     >               cdomvxy(ip,1,ivac),
     >               nmzxyd,ng2,kimax2,1,
     >               igfft2,pgfft2,nstr2)
            DO i=0,9*k1d*k2d-1 
              chdens= (af2(i,1)+af2(i,2))/2.  
              magmom= mx(i)**2 + my(i)**2 +
     &                ((af2(i,1)-af2(i,2))/2.)**2 
              magmom= SQRT(magmom) 
              af2(i,1)= chdens + magmom 
              af2(i,2)= chdens - magmom
            END DO 

          END IF 
c
c         calculate the exchange-correlation potential in  real space
c
          CALL vxcall
     >               (6,icorr,krla,jspins,
     >                ifftd2,nt,af2,
     <                vxc)   

          DO 80 js = 1,jspins
c
c            ----> 2-d back fft to g space
c
             CALL dset(ifftd2,0.0,bf2,1)
             CALL fft2d(
     >                 k1d,k2d,n2d,
     <                 vxc(0,js),bf2,
     >                 fgz,rhti,fgxy,
     >                 1,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)
c
c            ----> and add vxc to coulomb potential
c            the G||.eq.zero component is added to vz
c
             vz(ip,ivac,js) = fgz + vz(ip,ivac,js)
c
c            the G||.ne.zero components are added to vxy
c
             DO irec2 = 1,ng2-1
                vxy(ip,irec2,ivac,js) = vxy(ip,irec2,ivac,js) +
     +                                              fgxy(irec2)
             ENDDO   
  80      CONTINUE
c
ci        calculate the exchange-correlation energy density in  real space
c
          IF (total) THEN   
            CALL excall
     >                 (6,icorr,krla,jspins,
     >                  ifftd2,nt,af2,
     <                  exc)   
c
c     ----> 2-d back fft to g space
c
             CALL dset(ifftd2,0.0,bf2,1)
             CALL fft2d(
     >                 k1d,k2d,n2d,
     <                 exc,bf2,
     >                 excz(ip,ivac),rhti,excxy(ip,1,ivac),
     >                 nmzxyd,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)

           ENDIF

  110    CONTINUE
c
c        calculate vxc for z now beyond warping region 
c
         nmzdiff = nmz - nmzxy

         DO k=1,nmzdiff

           DO js=1,jspins
             af2(k-1,js) = rht(nmzxy+k,ivac,js)
           ENDDO

           IF (l_noco) THEN

             mx(0)= REAL(cdomvz(nmzxy+k,ivac))
             my(0)= AIMAG(cdomvz(nmzxy+k,ivac))
             chdens= (af2(k-1,1)+af2(k-1,2))/2.
             magmom= mx(0)**2 + my(0)**2 +
     &               ((af2(k-1,1)-af2(k-1,2))/2.)**2
             magmom= SQRT(magmom)
             af2(k-1,1)= chdens + magmom
             af2(k-1,2)= chdens - magmom

           END IF 

         ENDDO

         CALL vxcall
     >              (6,icorr,krla,jspins,
     >               nmzd,nmzdiff,af2,
     <               vxcz)
c+gu
         DO  js=1,jspins
           DO k=nmzxy+1,nmz
             vz(k,ivac,js) = vz(k,ivac,js) + vxcz(k-nmzxy,js)
           ENDDO
         ENDDO
c
         WRITE (6,FMT=8020) ivac, (vz(nmz,ivac,js),js=1,jspins)
         WRITE (16,FMT=8020) ivac, (vz(nmz,ivac,js),js=1,jspins)
 8020    FORMAT (/,5x,'vacuum zero for vacuum',i3,' = ',2f10.5)
c
c        calculate the ex.-corr. energy density now beyond warping region
c
         IF (total) THEN
            CALL excall
     >                   (6,icorr,krla,jspins,
     >                    nmzd,nmzdiff,af2,
     <                    excz(nmzxy+1,ivac)) 
         END IF
  150 CONTINUE

      IF (l_noco) THEN 
        DEALLOCATE (mx,my)
      END IF 

      END SUBROUTINE vvacxc
      END MODULE m_vvacxc
