      MODULE m_pldngen
c**********************************************************************
c     This subroutine generates the charge and magetization densities
c     (mx,my,mz) and writes them to the files cdn, mdnx, mdny, mdnz.
c     These files are needed to generate plots of the density.
c
c    i) The components of the hermitian density matrix (rho_11, rho_22,
c     rho_21) are reloaded from the file rhomat_inp.
c    ii) The density matrix in fouriertransformed to real space.
c    iii) The charge and magnetization density (n, mx, my, mz) are
c     calculated on the real space mesh.
c    iv) n, mx, my, and mz are Fouriertransformed and stored in terms
c     of stars.
c
c     Philipp Kurz 99/10/29
c**********************************************************************
      CONTAINS
      SUBROUTINE pldngen(
     >                   jspd,nop,k1d,k2d,k3d,n2d,n3d,
     >                   ntypd,ntypsd,nlhd,jmtd,natd,nmzd,nmzxyd,
     >                   nrhomfile,neq,invs,invs2,rmt,rmsh,zatom,namat,
     >                   dx,delz,z1,sk3,tau,bmat,ig,ig2,symor,vol,taual,
     >                   jspins,film,nvac,volmts,volint,area,
     >                   ntyp,jri,ntypsy,nlh,alph,beta,
     >                   mrot,ng2,ng3,kv3,igfft2,invtab,
     >                   igfft,pgfft2,pgfft,nstr2,nstr,
     >                   kimax,kimax2,nmz,nmzxy,sigma,odi,slice)

c******** ABBREVIATIONS ***********************************************
c     ifft3    : size of the 3d real space mesh
c     ifft2    : size of the 2d real space mesh
c     rpw      : first diagonal components of the interstitial density
c                matrix
c                later charge and mag. density (n, mx, my, mz)
c                all stored in terms of 3d-stars
c     ris      : first componets of the density matrix
c                later charge and mag. density (n, mx, my, mz)
c                all stored on real space mesh
c**********************************************************************

      USE m_loddop
      USE m_wrtdop
      USE m_qfix
      USE m_fft2d
      USE m_fft3d
      USE m_set,      ONLY : dset
      USE m_od_types, ONLY : od_inp
      USE m_rotdenmat 

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,nop,k1d,k2d,k3d,n2d,n3d,jmtd
      INTEGER, INTENT (IN) :: nrhomfile,ntypd,natd,ntypsd,nlhd
      INTEGER, INTENT (IN) :: jspins,nvac,ntyp,nmz,nmzxy,nmzd,nmzxyd
      INTEGER, INTENT (IN) :: ng2,ng3,kimax,kimax2
      REAL,    INTENT (IN) :: vol,volint,area,sigma,delz,z1
      LOGICAL, INTENT (IN) :: film,invs,invs2,symor,slice 
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),invtab(nop)
      INTEGER, INTENT (IN) :: nlh(ntypsd),neq(ntypd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: nstr2(n2d),nstr(n3d),kv3(3,n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL,    INTENT (IN) :: rmt(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),sk3(n3d),dx(ntypd),zatom(ntypd)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd)
      CHARACTER*2,INTENT (IN):: namat(0:103)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER iden,ivac,ifft2,ifft3
      INTEGER imz,ityp,iri,ilh,imesh,lh,iq2,iq3,iter
      REAL cdnup,cdndown,chden,mgden,theta,phi,zero,rho_11,rziw
      REAL rho_22,rho_21r,rho_21i,rhotot,mx,my,mz,fix,vz_r,vz_i
      COMPLEX czero
      CHARACTER*8 dop,iop,name(10)
C     ..
C     .. Local Arrays ..
c---> off-diagonal part of the density matrix
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      REAL,    ALLOCATABLE :: rht(:,:,:),rho(:,:,:,:)
      REAL,    ALLOCATABLE :: rvacxy(:,:,:,:),ris(:,:),fftwork(:)

c---> for testing: output of offdiag. output density matrix. to plot the
c---> offdiag. part of the output density matrix, that part has to be
c---> written the file rhomt21 in cdnmt.
      LOGICAL :: l_fmpl2
      REAL    :: cdn11, cdn22  
      COMPLEX :: cdn21 
      COMPLEX, ALLOCATABLE :: rho21(:,:,:)
c---> end of test part
c
      zero = 0.0 ; czero = cmplx(0.0,0.0)
      ifft3 = 27*k1d*k2d*k3d
      ifft2 = 9*k1d*k2d

      ALLOCATE (qpw(n3d,4),rhtxy(nmzxyd,n2d-1,2,4),
     +          cdom(n3d),cdomvz(nmzd,2),cdomvxy(nmzxyd,n2d-1,2),
     +     ris(0:27*k1d*k2d*k3d-1,4),fftwork(0:27*k1d*k2d*k3d-1),
     +     rvacxy(0:9*k1d*k2d-1,nmzxyd,2,4),
     +     rho(jmtd,0:nlhd,ntypd,4),rht(nmzd,2,4) )
!
!---> initialize arrays for the density matrix
!
      rho(:,:,:,:) = zero ; qpw(:,:) = czero ; cdom(:) = czero
      IF (film) THEN
        cdomvz(:,:) = czero ;    rhtxy(:,:,:,:) = czero
        cdomvxy(:,:,:) = czero ; rht(:,:,:) = zero
      ENDIF

      IF (jspins .NE. 2) THEN
         WRITE (6,*) 'This is the non-collinear version of the flapw-'
         WRITE (6,*) 'program. It can only perform spin-polarized'
         WRITE (6,*) 'calculations.'
         STOP 'jspins not equal 2'
      ENDIF

c---> reload the density matrix from file rhomat_inp
      IF ( .not. slice ) THEN  
        OPEN (nrhomfile,FILE='rhomat_inp',FORM='unformatted',
     +        STATUS='unknown')
      ELSE
        OPEN (nrhomfile,FILE='cdn_slice',FORM='unformatted',
     +        STATUS='unknown')
      ENDIF 
c---> first the diagonal elements of the density matrix
      CALL loddop(jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,ng3,ng2,nvac,ntyp,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nrhomfile,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c---> and then the off-diagonal part
      READ (nrhomfile,END=100,ERR=50) (cdom(iq3),iq3=1,ng3)
      IF (film) THEN
         READ (nrhomfile,END=75,ERR=50) ((cdomvz(imz,ivac),imz=1,nmz)
     +                              ,ivac=1,nvac)
         READ (nrhomfile,END=75,ERR=50) (((cdomvxy(imz,iq2-1,ivac)
     +                       ,imz=1,nmzxy),iq2=2,ng2),ivac=1,nvac)
      ENDIF
      GOTO 150
 50   WRITE(6,*)'rhodirgen: ERROR: Problems while reading density'
      WRITE(6,*)'matrix from file rhomat_inp.'
      STOP 'rhomatdir: ERROR while reading file rhomat_inp'
 75   WRITE(6,*)'rhomatdir: ERROR: reached end of file rhomat_inp'
      WRITE(6,*)'while reading the vacuum part of the off-diagonal'
      WRITE(6,*)'element of the desity matrix.'
      STOP 'rhomatdir: ERROR while reading file rhomat_inp'
 100  WRITE(6,*)'rhodirgen: WARNING: The file rhomat_inp does not'
      WRITE(6,*)'contain off-diagonal part of the density matrix.'
      WRITE(6,*)'Assuming collinear magnetization.'
 150  CLOSE (nrhomfile)
      IF (.not. slice) THEN 
        CALL qfix(
     >          k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >          nlhd,nmzd,nmz,jspins,film,nvac,area,ng3,nmzxy,n2d,
     >          ntyp,neq,volmts,taual,z1,vol,volint,ng2,invtab,
     >          symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >          nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     X          qpw,rhtxy,rho,rht,odi,
     <          fix)
      ENDIF 

c---> for testing: read offdiag. output density matrix
      INQUIRE (file= 'rhomt21', exist= l_fmpl2)
      IF (l_fmpl2) THEN
        ALLOCATE( rho21(jmtd,0:nlhd,ntypd) )
        OPEN (26,file='rhomt21',form='unformatted',status='unknown')
        READ (26) rho21
        CLOSE (26)
      ENDIF 
c---> end of test output

c---> calculate the charge and magnetization density in the muffin tins
      DO ityp = 1,ntyp
         DO ilh = 0,nlh(ntypsy(ityp))
            DO iri = 1,jri(ntyp)
               IF (.not. l_fmpl2) THEN 
                 cdnup   = rho(iri,ilh,ityp,1)
                 cdndown = rho(iri,ilh,ityp,2)
                 theta = beta(ityp)
                 phi   = alph(ityp)
                 chden  = cdnup + cdndown
                 mgden  = cdnup - cdndown
                 rho(iri,ilh,ityp,1) = chden
                 rho(iri,ilh,ityp,2) = mgden*cos(phi)*sin(theta)
                 rho(iri,ilh,ityp,3) = mgden*sin(phi)*sin(theta)
                 rho(iri,ilh,ityp,4) = mgden*cos(theta)
               ELSE 
c--->            for testing: output of offdiag. output density matrix
                 cdn11 = rho(iri,ilh,ityp,1)
                 cdn22 = rho(iri,ilh,ityp,2)
                 cdn21 = rho21(iri,ilh,ityp)
                 CALL rot_den_mat(alph(ityp),beta(ityp),
     X                            cdn11,cdn22,cdn21)
                 rho(iri,ilh,ityp,1) = cdn11 + cdn22
                 rho(iri,ilh,ityp,2) = 2*real(cdn21)
                 rho(iri,ilh,ityp,3) = 2*aimag(cdn21)
                 rho(iri,ilh,ityp,4) = cdn11 - cdn22
c--->            end of test part
               ENDIF 
            ENDDO
         ENDDO
      ENDDO

      IF (l_fmpl2) THEN
        DEALLOCATE( rho21 )
      ENDIF 

c---> fouriertransform the diagonal part of the density matrix
c---> in the interstitial, qpw, to real space (ris)
      DO iden = 1,2
         CALL fft3d(ris(0,iden),fftwork,qpw(1,iden),
     >               k1d,k2d,k3d,
     >               ng3,kimax,+1,
     >               igfft,pgfft,nstr)
      ENDDO
c---> fouriertransform the off-diagonal part of the density matrix
      CALL fft3d(ris(0,3),ris(0,4),cdom(1),
     >               k1d,k2d,k3d,
     >               ng3,kimax,+1,
     >               igfft,pgfft,nstr)

c---> calculate the charge and magnetization density on the
c---> real space mesh
      DO imesh = 0,ifft3-1
         rho_11  = ris(imesh,1)
         rho_22  = ris(imesh,2)
         rho_21r = ris(imesh,3)
         rho_21i = ris(imesh,4)
         rhotot  = rho_11 + rho_22
         mx      =  2*rho_21r
         my      = -2*rho_21i
         mz      = (rho_11-rho_22)
         
         ris(imesh,1) = rhotot
         ris(imesh,2) = mx
         ris(imesh,3) = my
         ris(imesh,4) = mz
      ENDDO

c---> Fouriertransform the density matrix back to reciprocal space
      DO iden = 1,4
         CALL dset(ifft3,zero,fftwork,1)
         CALL fft3d(ris(0,iden),fftwork,qpw(1,iden),
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
                  CALL fft2d(
     >                 k1d,k2d,n2d,
     =                 rvacxy(0,imz,ivac,iden),fftwork,
     =                 rht(imz,ivac,iden),rziw,rhtxy(imz,1,ivac,iden),
     >                 nmzxyd,ng2,kimax2,1,
     >                 igfft2,pgfft2,nstr2)
               ENDDO
            ENDDO
         ENDDO
c--->    fouriertransform the off-diagonal part of the density matrix
         DO ivac = 1,nvac
            DO imz = 1,nmzxyd
               rziw = 0.0
               vz_r = real(cdomvz(imz,ivac))
               vz_i = aimag(cdomvz(imz,ivac))
               CALL fft2d(
     >              k1d,k2d,n2d,
     =              rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),
     =              vz_r,vz_i,
     =              cdomvxy(imz,1,ivac),
     >              nmzxyd,ng2,kimax2,1,
     >              igfft2,pgfft2,nstr2)
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
                  rhotot  = rho_11 + rho_22
                  mx      =  2*rho_21r
                  my      = -2*rho_21i
                  mz      = (rho_11-rho_22)
                  
                  rvacxy(imesh,imz,ivac,1) = rhotot
                  rvacxy(imesh,imz,ivac,2) = mx
                  rvacxy(imesh,imz,ivac,3) = my
                  rvacxy(imesh,imz,ivac,4) = mz
               ENDDO
            ENDDO
            DO imz = nmzxyd+1,nmzd
               rho_11  = rht(imz,ivac,1)
               rho_22  = rht(imz,ivac,2)
               rho_21r = real(cdomvz(imz,ivac))
               rho_21i = aimag(cdomvz(imz,ivac))
               rhotot  = rho_11 + rho_22
               mx      =  2*rho_21r
               my      = -2*rho_21i
               mz      = (rho_11-rho_22)

               rht(imz,ivac,1) = rhotot
               rht(imz,ivac,2) = mx
               rht(imz,ivac,3) = my
               rht(imz,ivac,4) = mz
            ENDDO
         ENDDO
c--->    Fouriertransform the matrix potential back to reciprocal space
         DO iden = 1,4
            DO ivac = 1,nvac
               DO imz = 1,nmzxyd
                  CALL dset(ifft3,zero,fftwork,1)
                  CALL fft2d(
     >                 k1d,k2d,n2d,
     =                 rvacxy(0,imz,ivac,iden),fftwork,
     =                 rht(imz,ivac,iden),rziw,rhtxy(imz,1,ivac,iden),
     >                 nmzxyd,ng2,kimax2,-1,
     >                 igfft2,pgfft2,nstr2)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

c---> save charge density to file cdn
      OPEN (72,FILE='cdn',FORM='unformatted',STATUS='unknown')
      CALL wrtdop(
     >            1,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            1,ng3,ng2,nmzxy,nmz,nvac,ntyp,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,
     >            72,iop,dop,iter,rho(1,0,1,1),qpw(1,1),rht(1,1,1),
     >            rhtxy(1,1,1,1),name)
      CLOSE (72)

c---> save mx to file mdnx
      OPEN (72,FILE='mdnx',FORM='unformatted',STATUS='unknown')
      CALL wrtdop(
     >            1,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            1,ng3,ng2,nmzxy,nmz,nvac,ntyp,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,
     >            72,iop,dop,iter,rho(1,0,1,2),qpw(1,2),rht(1,1,2),
     >            rhtxy(1,1,1,2),name)
      CLOSE (72)

c---> save my to file mdny
      OPEN (72,FILE='mdny',FORM='unformatted',STATUS='unknown')
      CALL wrtdop(
     >            1,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            1,ng3,ng2,nmzxy,nmz,nvac,ntyp,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,
     >            72,iop,dop,iter,rho(1,0,1,3),qpw(1,3),rht(1,1,3),
     >            rhtxy(1,1,1,3),name)
      CLOSE (72)

c---> save mz to file mdnz
      OPEN (72,FILE='mdnz',FORM='unformatted',STATUS='unknown')
      CALL wrtdop(
     >            1,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            1,ng3,ng2,nmzxy,nmz,nvac,ntyp,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,
     >            72,iop,dop,iter,rho(1,0,1,4),qpw(1,4),rht(1,1,4),
     >            rhtxy(1,1,1,4),name)
      CLOSE (72)

      DEALLOCATE (qpw,rhtxy,cdom,cdomvz,cdomvxy,
     +            ris,fftwork,rvacxy,rho,rht)

      RETURN
      END SUBROUTINE pldngen
      END MODULE m_pldngen
