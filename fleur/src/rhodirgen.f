      MODULE m_rhodirgen
c**********************************************************************
c     This subroutine calculates the spin-up and -down density, which
c     are needed to calculate the potential and writes them to the file
c     cdn. The local angle of the magnetization is kept in real space
c     and written to the file dirofmag. This is done in four steps.
c
c    i) The components of the hermitian density matrix (rho_11, rho_22,
c     rho_21) are reloaded from the file rhomat_inp.
c    ii) The density matrix in fouriertransformed to real space.
c    iii) The spin-up and -down densities and the local angle of the
c     magnetization are calculated on the real space mesh.    
c    iv) The spin-up and -down densities are Fouriertransformed, stored
c     in terms of stars and written to the file cdn. The local angle of
c     magnetization is kept on the real space mesh and written to the
c     file dirofmag.
c
c     Philipp Kurz 99/11/01
c**********************************************************************
      CONTAINS
      SUBROUTINE rhodirgen(
     >                     jspd,nop,k1d,k2d,k3d,n2d,n3d,
     >                     ntypd,ntypsd,nlhd,jmtd,nmzd,nmzxyd,natd,
     >                     nrhomfile,ndomfile,neq,invs,invs2,rmt,rmsh,
     >                     zatom,namat,dx,delz,z1,sk3,tau,bmat,ig,ig2,
     >                     symor,vol,taual,jspins,film,nvac,volmts,
     >                     volint,area,ntyp,jri,ntypsy,nlh,
     >                     mrot,ng2,ng3,kv3,igfft2,invtab,
     >                     igfft,pgfft2,pgfft,nstr2,nstr,
     >                     kimax,kimax2,nmz,nmzxy,sigma,odi,odl)

c******** ABBREVIATIONS ***********************************************
c     ifft3    : size of the 3d real space mesh
c     ifft2    : size of the 2d real space mesh
c     rpw      : diagonal components of the density matrix (rho_11 ,
c                rho_22)
c                later interstitial spin-up and -down density
c                all stored in terms of 3d-stars
c     ris      : first components of the density matrix
c                later interstitial spin-up and -down density and
c                direction of magnetic field (theta and phi)
c                all stored on real space mesh
c**********************************************************************

      USE m_constants, ONLY : pimach
      USE m_loddop
      USE m_wrtdop
      USE m_qfix
      USE m_fft2d
      USE m_fft3d
      USE m_set, ONLY : dset
      USE m_od_types

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,nop,k1d,k2d,k3d,n2d,n3d,jmtd
      INTEGER, INTENT (IN) :: nrhomfile,ndomfile,ntypd,natd,ntypsd,nlhd
      INTEGER, INTENT (IN) :: jspins,nvac,ntyp,nmz,nmzxy,nmzd,nmzxyd
      INTEGER, INTENT (IN) :: ng2,ng3,kimax,kimax2
      REAL,    INTENT (IN) :: vol,volint,area,sigma,delz,z1
      LOGICAL, INTENT (IN) :: film,invs,invs2,symor
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: nlh(ntypsd),neq(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: nstr2(n2d),nstr(n3d),kv3(3,n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL,    INTENT (IN) :: rmt(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),sk3(n3d),dx(ntypd),zatom(ntypd)
      CHARACTER*2,INTENT (IN):: namat(0:103)
C     ..
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_lda), INTENT (IN) :: odl
c+odim
C     .. Local Scalars ..
      INTEGER iden,jspin,ivac,ifft2,ifft3
      INTEGER imz,ityp,iri,ilh,imesh,iq2,iq3,iter
      REAL theta,phi,zero,rho_11,rho_22,rho_21r,rho_21i,rhotot,magmom
      REAL rho_up,rho_down,mx,my,mz,eps,pi,fix,vz_r,vz_i,rziw
      COMPLEX czero
      CHARACTER*8 dop,iop,name(10)
C     ..
C     .. Local Arrays ..
c---> off-diagonal part of the density matrix
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      REAL,    ALLOCATABLE :: rht(:,:,:),rho(:,:,:,:),rz(:,:,:)
      REAL,    ALLOCATABLE :: rvacxy(:,:,:,:),ris(:,:),fftwork(:)
C     ..
      zero = 0.0 ; czero = cmplx(0.0,0.0) 
      eps = 1.0e-20
      pi = pimach()

      ALLOCATE (qpw(n3d,jspd),rhtxy(nmzxyd,odi%n2d-1,2,jspd),
     +          cdom(n3d),cdomvz(nmzd,2),cdomvxy(nmzxyd,odi%n2d-1,2),
     +     ris(0:27*k1d*k2d*k3d-1,4),fftwork(0:27*k1d*k2d*k3d-1),
     +     rz(nmzd,2,2),
     +     rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd) )
!
!---> initialize arrays for the density matrix
!
      rho(:,:,:,:) = zero ; qpw(:,:) = czero ; cdom(:) = czero
      IF (film) THEN
        rht(:,:,:) = zero   ; rz(:,:,:) = zero
        cdomvz(:,:) = czero ; rhtxy(:,:,:,:) = czero
        cdomvxy(:,:,:) = czero
      ENDIF

      ifft3 = 27*k1d*k2d*k3d
      ifft2 = 9*k1d*k2d
      if (odi%d1) ifft2 = 9*k3d*odi%M
      if (film) allocate(rvacxy(0:ifft2-1,nmzxyd,2,4))

      IF (jspins .NE. 2) THEN
         WRITE (6,*) 'This is the non-collinear version of the flapw-'
         WRITE (6,*) 'program. It can only perform spin-polarized'
         WRITE (6,*) 'calculations.'
         STOP 'jspins not equal 2'
      ENDIF

c---> reload the density matrix from file rhomat_inp
      OPEN (nrhomfile,FILE='rhomat_inp',FORM='unformatted',
     +      STATUS='unknown')
c---> first the diagonal elements of the density matrix
      CALL loddop(jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,ng3,odi%nq2,nvac,ntyp,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nrhomfile,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c---> and then the off-diagonal part
      READ (nrhomfile,END=100,ERR=50) (cdom(iq3),iq3=1,ng3)
      IF (film) THEN
         READ (nrhomfile,END=75,ERR=50) ((cdomvz(imz,ivac),imz=1,nmz)
     +                              ,ivac=1,nvac)
         READ (nrhomfile,END=75,ERR=50) (((cdomvxy(imz,iq2-1,ivac)
     +                       ,imz=1,nmzxy),iq2=2,odi%nq2),ivac=1,nvac)
      ENDIF
      GOTO 150
 50   WRITE(6,*)'rhodirgen: ERROR: Problems while reading density'
      WRITE(6,*)'matrix from file rhomat_inp.'
      STOP 'rhodirgen: ERROR while reading file rhomat_inp'
 75   WRITE(6,*)'rhodirgen: ERROR: reached end of file rhomat_inp'
      WRITE(6,*)'while reading the vacuum part of the off-diagonal'
      WRITE(6,*)'element of the desity matrix.'
      STOP 'rhodirgen: ERROR while reading file rhomat_inp'
 100  WRITE(6,*)'rhodirgen: WARNING: The file rhomat_inp does not'
      WRITE(6,*)'contain off-diagonal part of the density matrix.'
      WRITE(6,*)'Assuming collinear magnetization.'
 150  CLOSE (nrhomfile)
      CALL qfix(
     >          k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >          nlhd,nmzd,nmz,jspins,film,nvac,area,ng3,nmzxy,n2d,
     >          ntyp,neq,volmts,taual,z1,vol,volint,ng2,invtab,
     >          symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >          nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     X          qpw,rhtxy,rho,rht,odi,
     <          fix)

c---> fouriertransform the diagonal part of the density matrix
c---> in the interstitial, qpw, to real space (ris)
      DO iden = 1,2
         CALL fft3d(
     <               ris(0,iden),fftwork,
     >               qpw(1,iden),
     >               k1d,k2d,k3d,
     >               ng3,kimax,+1,
     >               igfft,pgfft,nstr)
      ENDDO
c---> fouriertransform the off-diagonal part of the density matrix
      CALL fft3d(
     <           ris(0,3),ris(0,4),
     >           cdom(1),
     >           k1d,k2d,k3d,
     >           ng3,kimax,+1,
     >           igfft,pgfft,nstr)

ctest
c      DO iden=1,4
c         write(*,*)'iden=',iden
c         write(*,8500)(ris(imesh,iden),imesh=0,ifft3-1)
c      enddo
ctest
c---> calculate the charge and magnetization density on the
c---> real space mesh
      DO imesh = 0,ifft3-1
         rho_11  = ris(imesh,1)
         rho_22  = ris(imesh,2)
         rho_21r = ris(imesh,3)
         rho_21i = ris(imesh,4)
         mx      =  2*rho_21r
         my      = -2*rho_21i
         mz      = (rho_11-rho_22)
         magmom  = sqrt(mx**2 + my**2 + mz**2)
         rhotot  = rho_11 + rho_22
         rho_up  = (rhotot + magmom)/2
         rho_down= (rhotot - magmom)/2

         IF (abs(mz) .LE. eps) THEN
            theta = pi/2
         ELSEIF (mz .GE. 0.0) THEN
            theta = atan(sqrt(mx**2 + my**2)/mz)
         ELSE
            theta = atan(sqrt(mx**2 + my**2)/mz) + pi
         ENDIF
         
         IF (abs(mx) .LE. eps) THEN
            IF (abs(my) .LE. eps) THEN
               phi = 0.0
            ELSEIF (my .GE. 0.0) THEN
               phi = pi/2
            ELSE
               phi = -pi/2
            ENDIF
         ELSEIF (mx .GE. 0.0) THEN
            phi = atan(my/mx)
         ELSE
            IF (my .GE. 0.0) THEN
               phi = atan(my/mx) + pi
            ELSE
               phi = atan(my/mx) - pi
            ENDIF
         ENDIF

c         write(36,'(i4,2f12.6)') mod(imesh,33),rho_11,rho_22
         ris(imesh,1) = rho_up
         ris(imesh,2) = rho_down
         ris(imesh,3) = theta
         ris(imesh,4) = phi
      ENDDO
ctest
c      DO iden=1,4
c         write(*,*)'iden=',iden
c         write(*,8500)(ris(imesh,iden),imesh=0,ifft3-1)
c 8500    format(10e13.5)
c      enddo
ctest
c---> Fouriertransform the density matrix back to reciprocal space
      DO jspin = 1,jspins
         CALL dset(ifft3,zero,fftwork,1)
         CALL fft3d(
     >               ris(0,jspin),fftwork,
     <               qpw(1,jspin),
     >               k1d,k2d,k3d,
     >               ng3,kimax,-1,
     >               igfft,pgfft,nstr)
      ENDDO

c---> fouriertransform the diagonal part of the density matrix
c---> in the vacuum, rz & rxy, to real space (rvacxy)
      IF (film) THEN
         DO iden = 1,2
            DO ivac = 1,nvac
               DO imz = 1,nmzxyd
                  rziw = 0.0
                  IF (odi%d1) THEN
                   CALL fft2d(
     >                 odi%k3,odi%M,odi%n2d,
     =                 rvacxy(0,imz,ivac,iden),fftwork,
     =                 rht(imz,ivac,iden),rziw,rhtxy(imz,1,ivac,iden),
     >                 nmzxyd,odi%nq2,odi%kimax2,1,
     >                 odl%igf,odl%pgf,odi%nst2)
                  ELSE
                   CALL fft2d(
     >                 k1d,k2d,n2d,
     =                 rvacxy(0,imz,ivac,iden),fftwork,
     =                 rht(imz,ivac,iden),rziw,rhtxy(imz,1,ivac,iden),
     >                 nmzxyd,ng2,kimax2,1,
     >                 igfft2,pgfft2,nstr2)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
c--->    fouriertransform the off-diagonal part of the density matrix
         DO ivac = 1,nvac
            DO imz = 1,nmzxyd
               rziw = 0.0
               vz_r = real(cdomvz(imz,ivac))
               vz_i = aimag(cdomvz(imz,ivac))
               IF (odi%d1) THEN
                CALL fft2d(
     >              odi%k3,odi%M,odi%n2d,
     =              rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),
     =              vz_r,vz_i,
     =              cdomvxy(imz,1,ivac),
     >              nmzxyd,odi%nq2,odi%kimax2,1,
     >              odl%igf,odl%pgf,odi%nst2)
               ELSE
                CALL fft2d(
     >              k1d,k2d,n2d,
     =              rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),
     =              vz_r,vz_i,
     =              cdomvxy(imz,1,ivac),
     >              nmzxyd,ng2,kimax2,1,
     >              igfft2,pgfft2,nstr2)
               ENDIF
            ENDDO
         ENDDO

c--->    calculate the four components of the matrix potential on
c--->    real space mesh
         DO ivac = 1,nvac
            DO imz = 1,nmzxyd
               DO imesh = 0,ifft2-1
                  rho_11  = rvacxy(imesh,imz,ivac,1)
                  rho_22  = rvacxy(imesh,imz,ivac,2)
                  rho_21r = rvacxy(imesh,imz,ivac,3)
                  rho_21i = rvacxy(imesh,imz,ivac,4)
                  mx      =  2*rho_21r
                  my      = -2*rho_21i
                  mz      = (rho_11-rho_22)
                  magmom  = sqrt(mx**2 + my**2 + mz**2)
                  rhotot  = rho_11 + rho_22
                  rho_up  = (rhotot + magmom)/2
                  rho_down= (rhotot - magmom)/2
                  
                  IF (abs(mz) .LE. eps) THEN
                     theta = pi/2
                  ELSEIF (mz .GE. 0.0) THEN
                     theta = atan(sqrt(mx**2 + my**2)/mz)
                  ELSE
                     theta = atan(sqrt(mx**2 + my**2)/mz) + pi
                  ENDIF

                  IF (abs(mx) .LE. eps) THEN
                     IF (abs(my) .LE. eps) THEN
                        phi = 0.0
                     ELSEIF (my .GE. 0.0) THEN
                        phi = pi/2
                     ELSE
                        phi = -pi/2
                     ENDIF
                  ELSEIF (mx .GE. 0.0) THEN
                     phi = atan(my/mx)
                  ELSE
                     IF (my .GE. 0.0) THEN
                        phi = atan(my/mx) + pi
                     ELSE
                        phi = atan(my/mx) - pi
                     ENDIF
                  ENDIF

                  rvacxy(imesh,imz,ivac,1) = rho_up
                  rvacxy(imesh,imz,ivac,2) = rho_down
                  rvacxy(imesh,imz,ivac,3) = theta
                  rvacxy(imesh,imz,ivac,4) = phi
               ENDDO
            ENDDO
            DO imz = nmzxyd+1,nmzd
               rho_11  = rht(imz,ivac,1)
               rho_22  = rht(imz,ivac,2)
               rho_21r = real(cdomvz(imz,ivac))
               rho_21i = aimag(cdomvz(imz,ivac))
               mx      =  2*rho_21r
               my      = -2*rho_21i
               mz      = (rho_11-rho_22)
               magmom  = sqrt(mx**2 + my**2 + mz**2)
               rhotot  = rho_11 + rho_22
               rho_up  = (rhotot + magmom)/2
               rho_down= (rhotot - magmom)/2

               IF (abs(mz) .LE. eps) THEN
                  theta = pi/2
               ELSEIF (mz .GE. 0.0) THEN
                  theta = atan(sqrt(mx**2 + my**2)/mz)
               ELSE
                  theta = atan(sqrt(mx**2 + my**2)/mz) + pi
               ENDIF
               
               IF (abs(mx) .LE. eps) THEN
                  IF (abs(my) .LE. eps) THEN
                     phi = 0.0
                  ELSEIF (my .GE. 0.0) THEN
                     phi = pi/2
                  ELSE
                     phi = -pi/2
                  ENDIF
               ELSEIF (mx .GE. 0.0) THEN
                  phi = atan(my/mx)
               ELSE
                  IF (my .GE. 0.0) THEN
                     phi = atan(my/mx) + pi
                  ELSE
                     phi = atan(my/mx) - pi
                  ENDIF
               ENDIF

               rht(imz,ivac,1) = rho_up
               rht(imz,ivac,2) = rho_down
               rz(imz,ivac,1) = theta
               rz(imz,ivac,2) = phi
            ENDDO
         ENDDO
c--->    Fouriertransform the matrix potential back to reciprocal space
         DO jspin = 1,jspins
            DO ivac = 1,nvac
               DO imz = 1,nmzxyd
                  CALL dset(ifft3,zero,fftwork,1)
                  IF (odi%d1) THEN
                   CALL fft2d(
     >                 odi%k3,odi%M,odi%n2d,
     =                 rvacxy(0,imz,ivac,jspin),fftwork,
     =                 rht(imz,ivac,jspin),rziw,rhtxy(imz,1,ivac,jspin),
     >                 nmzxyd,odi%nq2,odi%kimax2,-1,
     >                 odl%igf,odl%pgf,odi%nst2)
                  ELSE
                   CALL fft2d(
     >                 k1d,k2d,n2d,
     =                 rvacxy(0,imz,ivac,jspin),fftwork,
     =                 rht(imz,ivac,jspin),rziw,rhtxy(imz,1,ivac,jspin),
     >                 nmzxyd,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)
                  END IF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c---> ndomfile is the dirofmag-file
      OPEN (ndomfile,FILE='dirofmag',FORM='unformatted',
     +     STATUS='unknown')
      WRITE (ndomfile) (ris(imesh,3),imesh=0,ifft3-1)
      WRITE (ndomfile) (ris(imesh,4),imesh=0,ifft3-1)
      IF (film) THEN
         WRITE (ndomfile) ((rz(imz,ivac,1),imz=nmzxyd+1,nmzd),
     +        ivac=1,nvac)
         WRITE (ndomfile) ((rz(imz,ivac,2),imz=nmzxyd+1,nmzd),
     +        ivac=1,nvac)
         WRITE (ndomfile) (((rvacxy(imesh,imz,ivac,3),
     +        imesh=0,ifft2-1),imz=1,nmzxyd),ivac=1,nvac)
         WRITE (ndomfile) (((rvacxy(imesh,imz,ivac,4),
     +        imesh=0,ifft2-1),imz=1,nmzxyd),ivac=1,nvac)
      ENDIF
      CLOSE (ndomfile)

c---> write spin-up and -down density on file cdn
      OPEN (70,FILE='cdn',FORM='unformatted',STATUS='unknown')
      CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,ng3,odi%nq2,nmzxy,nmz,nvac,ntyp,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,
     >            70,iop,dop,iter,rho,qpw,rht,rhtxy,name)
      CLOSE (70)

      DEALLOCATE (qpw,rhtxy,cdom,cdomvz,cdomvxy,
     +            ris,fftwork,rz,rho,rht)
      if (film) deallocate(rvacxy)
      RETURN
      END SUBROUTINE rhodirgen
      END MODULE m_rhodirgen
