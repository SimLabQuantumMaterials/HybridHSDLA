      MODULE m_vmatgen
c**********************************************************************
c     This subroutine prepares the spin dependent 2x2 matrix potential
c     for the Hamiltonian setup. This is done in 4 steps.
c
c    i) The spin up and down potential and the direction of the
c     magentic field, theta and phi, are reloaded from files nrp,
c     dirofmag.
c    ii) The spin up and down potential is Fouriertransformed to real
c     space (theta and phi are stored in real space).
c    iii) The four components of the matrix potential are calculated on
c     the real space mesh.
c    iv) The matrix potential is Fouriertransformed, stored in terms of
c     stars and written to file potmat.
c
c     Philipp Kurz 99/11/01
c**********************************************************************
      CONTAINS
      SUBROUTINE vmatgen(
     >                   jspd,k1d,k2d,k3d,n2d,n3d,
     >                   ntypd,ntypsd,nlhd,natd,jmtd,nmzd,nmzxyd,
     >                   ntyp,ntypsy,nlh,jri,neq,invs,invs2,
     >                   nu,ndomfile,npotmatfile,jspins,film,nvac,
     >                   ng2,ng3,igfft2,
     >                   igfft,pgfft2,pgfft,nstr2,nstr,
     >                   kimax,kimax2,ufft,odi,odl)

c******** ABBREVIATIONS ***********************************************
c     ifft3    : size of the 3d real space mesh
c     ifft2    : size of the 2d real space mesh
c     vpw      : first interstitial spin up and down potential
c                later four components of matrix potential
c                all stored in terms of 3d-stars
c     vis      : first interstitial spin up and down potential and
c                direction of magnetic field (theta and phi)
c                later four components of matrix potential
c                all stored on real space mesh
c**********************************************************************

      USE m_loddop
      USE m_fft2d
      USE m_fft3d
      USE m_set, ONLY : dset
      USE m_od_types
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,k1d,k2d,k3d,n2d,n3d
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nlhd,natd,jmtd,nmzd,nmzxyd
      INTEGER, INTENT (IN) :: ntyp,ng2,ng3,kimax,kimax2
      INTEGER, INTENT (IN) :: nu,ndomfile,npotmatfile,jspins,nvac
      LOGICAL, INTENT (IN) :: film,invs,invs2

C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nstr2(n2d),nstr(n3d),neq(ntypd)
      INTEGER, INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL   , INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL   , INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL   , INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_lda), INTENT (IN) :: odl
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER imeshpt,ipot,jspin,ig2,ig3,ivac,ifft2,ifft3,imz,iter
      REAl    vup,vdown,veff,beff,theta,phi,zero,vziw
      LOGICAL l_domfexst
      CHARACTER*8 dop,iop,name(10)
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: vpw(:,:),vxy(:,:,:,:)
      REAL,    ALLOCATABLE :: vr(:,:,:,:),vz(:,:,:)
      REAL,    ALLOCATABLE :: vvacxy(:,:,:,:),vis(:,:),fftwork(:)

      zero = 0.0
      ifft3 = 27*k1d*k2d*k3d
      ifft2 = 9*k1d*k2d
      if (odi%d1) ifft2 = 9*k3d*odi%M
      if (film) allocate(vvacxy(0:ifft2-1,nmzxyd,2,4))

      IF (jspins .NE. 2) THEN
         WRITE (6,*) 'This is the non-collinear version of the flapw-'
         WRITE (6,*) 'program. It can only perform spin-polarized'
         WRITE (6,*) 'calculations.'
         STOP 'jspins not equal 2'
      ENDIF

      ALLOCATE ( vpw(n3d,3),vis(0:27*k1d*k2d*k3d-1,4),
     +           vxy(nmzxyd,odi%n2d-1,2,3),vr(jmtd,0:nlhd,ntypd,jspd),
     +           vz(nmzd,2,4),fftwork(0:27*k1d*k2d*k3d-1) )

c---> reload the spin up and down potential
c      OPEN (nu,file='pottot',form='unformatted',status='old')
      OPEN (nu,file='nrp',form='unformatted',status='old')
      CALL loddop(jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,ng3,odi%nq2,nvac,ntyp,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nu,natd,neq,
     <            iop,dop,iter,vr,vpw(1,1),vz(1,1,1),vxy(1,1,1,1),name)
      CLOSE(nu)

c---> check, whether the direction of magnetic field file exists
      INQUIRE (FILE='dirofmag',EXIST=l_domfexst)
      IF (l_domfexst) THEN
c--->    if it does, read the theta and phi values
         OPEN (ndomfile,FILE='dirofmag',FORM='unformatted',
     +        STATUS='unknown')
         READ (ndomfile) (vis(imeshpt,3),imeshpt=0,ifft3-1)
         READ (ndomfile) (vis(imeshpt,4),imeshpt=0,ifft3-1)
         IF (film) THEN
            READ (ndomfile) ((vz(imz,ivac,3),imz=nmzxyd+1,nmzd),
     +                       ivac=1,nvac)
            READ (ndomfile) ((vz(imz,ivac,4),imz=nmzxyd+1,nmzd),
     +                       ivac=1,nvac)
            READ (ndomfile) (((vvacxy(imeshpt,imz,ivac,3),
     +                   imeshpt=0,ifft2-1),imz=1,nmzxyd),ivac=1,nvac)
            READ (ndomfile) (((vvacxy(imeshpt,imz,ivac,4),
     +                   imeshpt=0,ifft2-1),imz=1,nmzxyd),ivac=1,nvac)
         ENDIF
         CLOSE (ndomfile)
      ELSE
c--->    if it doesn't, set all angles to zero
         CALL dset(ifft3,0.0,vis(0,3),1)
         CALL dset(ifft3,0.0,vis(0,4),1)
         IF (film) THEN
            DO ivac = 1,2
               DO imz = nmzxyd+1,nmzd
                  vz(imz,ivac,3) = 0.0
                  vz(imz,ivac,4) = 0.0
               ENDDO
               DO imz = 1,nmzxyd
                  DO imeshpt = 0,ifft2-1
                     vvacxy(imeshpt,imz,ivac,3) = 0.0
                     vvacxy(imeshpt,imz,ivac,4) = 0.0
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

c---> fouriertransform the spin up and down potential
c---> in the interstitial, vpw, to real space (vis)
      DO jspin = 1,jspins
         CALL fft3d(
     <               vis(0,jspin),fftwork,
     >               vpw(1,jspin),
     >               k1d,k2d,k3d,ng3,kimax,+1,
     >               igfft,pgfft,nstr)
      ENDDO

c---> calculate the four components of the matrix potential on
c---> real space mesh
      DO imeshpt = 0,ifft3-1
         vup   = vis(imeshpt,1)
         vdown = vis(imeshpt,2)
         theta = vis(imeshpt,3)
         phi   = vis(imeshpt,4)
c         write (35,'(i4,4f12.6)') mod(imeshpt,33),vup,vdown,theta,phi
c--->    at first determine the effective potential and magnetic field
         veff  = (vup + vdown)/2.0
         beff  = (vup - vdown)/2.0
c--->    now calculate the matrix potential, which is hermitian.
c--->    thus calculate the diagonal elements:
c--->    V_11
         vis(imeshpt,1) = veff + beff*cos(theta)
c--->    V_22
         vis(imeshpt,2) = veff - beff*cos(theta)
c--->    the real part of V_21
         vis(imeshpt,3) = beff*sin(theta)*cos(phi)
c--->    and the imaginary part of V_21
         vis(imeshpt,4) = beff*sin(theta)*sin(phi)
         DO ipot = 1,4
           vis(imeshpt,ipot) =  vis(imeshpt,ipot) * ufft(imeshpt)
         ENDDO
      ENDDO

c---> Fouriertransform the matrix potential back to reciprocal space
      DO ipot = 1,2
         CALL dset(ifft3,zero,fftwork,1)
         CALL fft3d(
     >               vis(0,ipot),fftwork,
     <               vpw(1,ipot),
     >               k1d,k2d,k3d,ng3,kimax,-1, 
     >               igfft,pgfft,nstr)
      ENDDO
      CALL fft3d(
     >           vis(0,3),vis(0,4),
     <           vpw(1,3),
     >           k1d,k2d,k3d,ng3,kimax,-1,
     >           igfft,pgfft,nstr)

      IF (film) THEN
c--->    fouriertransform the spin up and down potential
c--->    in the vacuum, vz & vxy, to real space (vvacxy)
         DO jspin = 1,jspins
            DO ivac = 1,nvac
               DO imz = 1,nmzxyd
                  vziw = 0.0
                  IF (odi%d1) THEN
                   CALL fft2d(
     >                 odi%k3,odi%M,odi%n2d,
     =                 vvacxy(0,imz,ivac,jspin),fftwork,
     =                 vz(imz,ivac,jspin),vziw,vxy(imz,1,ivac,jspin),
     >                 nmzxyd,odi%nq2,odi%kimax2,1,
     >                 odl%igf,odl%pgf,odi%nst2)
                  ELSE
                   CALL fft2d(
     >                 k1d,k2d,n2d,
     =                 vvacxy(0,imz,ivac,jspin),fftwork,
     =                 vz(imz,ivac,jspin),vziw,vxy(imz,1,ivac,jspin),
     >                 nmzxyd,ng2,kimax2,1,
     >                 igfft2,pgfft2,nstr2)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

c--->    calculate the four components of the matrix potential on
c--->    real space mesh
         DO ivac = 1,nvac
            DO imz = 1,nmzxyd
               DO imeshpt = 0,ifft2-1
                  vup   = vvacxy(imeshpt,imz,ivac,1)
                  vdown = vvacxy(imeshpt,imz,ivac,2)
                  theta = vvacxy(imeshpt,imz,ivac,3)
                  phi   = vvacxy(imeshpt,imz,ivac,4)
                  veff  = (vup + vdown)/2.0
                  beff  = (vup - vdown)/2.0
                  vvacxy(imeshpt,imz,ivac,1) = veff + beff*cos(theta)
                  vvacxy(imeshpt,imz,ivac,2) = veff - beff*cos(theta)
                  vvacxy(imeshpt,imz,ivac,3) = beff*sin(theta)*cos(phi)
                  vvacxy(imeshpt,imz,ivac,4) = beff*sin(theta)*sin(phi)
               ENDDO
            ENDDO
            DO imz = nmzxyd+1,nmzd
               vup   = vz(imz,ivac,1)
               vdown = vz(imz,ivac,2)
               theta = vz(imz,ivac,3)
               phi   = vz(imz,ivac,4)
               veff  = (vup + vdown)/2.0
               beff  = (vup - vdown)/2.0
               vz(imz,ivac,1) = veff + beff*cos(theta)
               vz(imz,ivac,2) = veff - beff*cos(theta)
               vz(imz,ivac,3) = beff*sin(theta)*cos(phi)
               vz(imz,ivac,4) = beff*sin(theta)*sin(phi)
            ENDDO
         ENDDO

c--->    Fouriertransform the matrix potential back to reciprocal space
         DO ipot = 1,2
            DO ivac = 1,nvac
               DO imz = 1,nmzxyd
                  CALL dset(ifft3,zero,fftwork,1)
                  IF (odi%d1) THEN
                   CALL fft2d(
     >                 odi%k3,odi%M,odi%n2d,
     =                 vvacxy(0,imz,ivac,ipot),fftwork,
     =                 vz(imz,ivac,ipot),vziw,vxy(imz,1,ivac,ipot),
     >                 nmzxyd,odi%nq2,odi%kimax2,-1,
     >                 odl%igf,odl%pgf,odi%nst2)
                  ELSE
                   CALL fft2d(
     >                 k1d,k2d,n2d,
     =                 vvacxy(0,imz,ivac,ipot),fftwork,
     =                 vz(imz,ivac,ipot),vziw,vxy(imz,1,ivac,ipot),
     >                 nmzxyd,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)
                  END IF
               ENDDO
            ENDDO
         ENDDO

         DO ivac = 1,nvac
            DO imz = 1,nmzxyd
               CALL dset(ifft3,zero,fftwork,1)
               IF (odi%d1) THEN
                CALL fft2d(
     >              odi%k3,odi%M,odi%n2d,
     =              vvacxy(0,imz,ivac,3),vvacxy(0,imz,ivac,4),
     =              vz(imz,ivac,3),vz(imz,ivac,4),vxy(imz,1,ivac,3),
     >              nmzxyd,odi%nq2,odi%kimax2,-1,
     >              odl%igf,odl%pgf,odi%nst2)
               ELSE
                CALL fft2d(
     >              k1d,k2d,n2d,
     =              vvacxy(0,imz,ivac,3),vvacxy(0,imz,ivac,4),
     =              vz(imz,ivac,3),vz(imz,ivac,4),vxy(imz,1,ivac,3),
     >              nmzxyd,ng2,kimax2,-1,
     >              igfft2,pgfft2,nstr2)
               END IF
            ENDDO
         ENDDO

      ENDIF
c
c---> save matrix potential to file potmat
c
      OPEN (npotmatfile,FILE='potmat',FORM='unformatted',
     +     STATUS='unknown')
      DO ipot = 1,3
        WRITE (npotmatfile) (vpw(ig3,ipot),ig3=1,ng3)
      ENDDO
      IF (film) THEN
         DO ivac = 1,nvac
            WRITE (npotmatfile)((vz(imz,ivac,ipot),imz=1,nmzd),ipot=1,4)
            DO ipot = 1,3
               WRITE (npotmatfile)((vxy(imz,ig2,ivac,ipot),
     +                      imz=1,nmzxyd),ig2=1,odi%nq2-1)
            ENDDO
         ENDDO
      ENDIF
 8000 FORMAT(6f16.10)
      CLOSE (npotmatfile)

      DEALLOCATE ( vpw,vis,vxy,vr,vz,fftwork)
      if (film) deallocate (vvacxy)
      RETURN
      END SUBROUTINE vmatgen
      END MODULE m_vmatgen
