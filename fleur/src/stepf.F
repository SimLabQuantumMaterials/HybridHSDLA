      MODULE m_stepf
      CONTAINS
      SUBROUTINE stepf(
     >                 k1d,k2d,k3d,n3d,natd,ntypd,odd,sk2,n2d,
     >                 nq3,ntype,film,neq,kv3,sk3,ig2,omtil,ig,
     >                 dvac,area,z1,vol,rmt,taual,bmat,volmts,
     <                 ufft,ustep)
c
c*********************************************************************
c     calculates the fourier components of the interstitial step
c     function for the reciprocal vectors of the star list.
c           m. weinert  1986
c*********************************************************************
c
c     also set up FFT of U(G) on a (-2G:+2G) grid for convolutions
c
c*********************************************************************

#include"cpp_double.h"
      USE m_constants, ONLY : pimach
      USE m_od_cylbes
      USE m_set, ONLY : dset
      USE m_od_types, ONLY : od_dim
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,natd,ntypd
      INTEGER, INTENT (IN) :: nq3,ntype,n2d
      LOGICAL, INTENT (IN) :: film
      REAL,    INTENT (IN) :: area,dvac,z1,vol,omtil
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),kv3(3,n3d),ig2(n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: sk3(n3d),rmt(ntypd),volmts(ntypd)
      REAL,    INTENT (IN) :: taual(3,natd),bmat(3,3),sk2(n2d)
      REAL,    INTENT (OUT) :: ufft(0:27*k1d*k2d*k3d-1)
      COMPLEX, INTENT (OUT) :: ustep(n3d)
c-odim
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX c_c,c_phs
      REAL c,dd,gs,th,inv_omtil,tpi,fpi,r_phs
      REAL g_rmt,g_sqr,help,g_abs,fp_omtil,r_c,gr,gx,gy
      INTEGER i,k,n,n3,na,nn,i1,i2,i3,ic,ifft2d,ifftd,kk
      INTEGER ic1,ic2,ic3,icc,im1,im2,im3
C     ..
C     .. Local Arrays ..
      COMPLEX sf(n3d)
      REAL g(3),gm(3),fJ
      REAL,    ALLOCATABLE :: bfft(:)
      INTEGER, ALLOCATABLE :: icm(:,:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,cos,exp,sin
C     ..
C     .. External Functions ..
      REAL CPP_BLAS_sdot
      EXTERNAL CPP_BLAS_sdot
C     ..
C     ..
c--->    if step function on unit14, then just read it in
c
      ifftd = 27*k1d*k2d*k3d
c
      OPEN (14,file='wkf2',form='unformatted',status='unknown')
      REWIND 14
      READ (14,end=10,err=10) n3,n
      IF (n3.NE.nq3) GO TO 10
      IF (n.NE.ifftd) GO TO 10
      READ (14) (ustep(i),i=1,nq3)
      READ (14) (ufft(i),i=0,ifftd-1)
      CLOSE (14)
      RETURN

   10 CONTINUE
      tpi = 2 * pimach()
      fpi = 2 * tpi
c
      IF (film) THEN
         dd = dvac*area/omtil
         IF (odd%d1) dd = vol/omtil
      ELSE
         dd = 1.0
      END IF
c--->    G=0 star
      c = 0.0
      DO 20 n = 1,ntype
         c = c + neq(n)*volmts(n)/omtil
   20 CONTINUE
      ustep(1) = cmplx(dd-c,0.0)
c--->    G(parallel)=0  (for film)
      IF (film .AND. .NOT.odd%d1) THEN
         DO 30 k = 2,nq3
            IF (ig2(k).EQ.1) THEN
               th = bmat(3,3)*kv3(3,k)*z1
               ustep(k) = cmplx(vol*sin(th)/th/omtil,0.0)
            ELSE
               ustep(k) = cmplx(0.0,0.0)
            END IF
   30    CONTINUE
c-odim
      ELSEIF (odd%d1) THEN
         DO 35 k = 2,nq3
            gr = 0.0
            IF (kv3(3,k).EQ.0) THEN
               kk = ig2(k)
               gr = sk2(kk)
               CALL od_cylbes(1,gr*z1,fJ)
               ustep(k) = cmplx(2.*dd*fJ/(gr*z1),0.)
            ELSE
               ustep(k) =cmplx(0.,0.)
            END IF
 35      CONTINUE
c+odim
      ELSE
         DO 40 k = 2,nq3
            ustep(k) = cmplx(0.0,0.0)
   40    CONTINUE
      END IF
c--->    sphere contributions
      na = 0
      DO 90 n = 1,ntype
         c = 3.*volmts(n)/omtil
c-->     structure factors: loop over equivalent atoms
         na = na + 1
         DO 50 k = 2,nq3
            th = -tpi* (kv3(1,k)*taual(1,na)+kv3(2,k)*taual(2,na)+
     +           kv3(3,k)*taual(3,na))
            sf(k) = cmplx(cos(th),sin(th))
   50    CONTINUE
         DO 70 nn = 2,neq(n)
            na = na + 1
            DO 60 k = 2,nq3
               th = -tpi* (kv3(1,k)*taual(1,na)+kv3(2,k)*taual(2,na)+
     +              kv3(3,k)*taual(3,na))
               sf(k) = sf(k) + cmplx(cos(th),sin(th))
   60       CONTINUE
   70    CONTINUE
c--->    update step function
         DO 80 k = 2,nq3
            gs = sk3(k)*rmt(n)
            ustep(k) = ustep(k) - (c* (sin(gs)/gs-cos(gs))/ (gs*gs))*
     +                 sf(k)
   80    CONTINUE
   90 CONTINUE
c
c --> set up stepfunction on fft-grid:
c
      ALLOCATE (  bfft(0:27*k1d*k2d*k3d-1) )
      im1=CEILING(1.5*k1d); im2=CEILING(1.5*k2d); im3=CEILING(1.5*k3d) 
      ALLOCATE ( icm(-im1:im1,-im2:im2,-im3:im3) )
      icm = 0
      ic=0
      inv_omtil=1.0/omtil
      fp_omtil=  -fpi*inv_omtil
      DO i3=0,3*k3d-1
        gm(3)=real(i3)
        IF ( gm(3) > 1.5*k3d ) gm(3)=gm(3)-3.0*k3d
        DO i2=0,3*k2d-1
          gm(2)=real(i2)
          IF ( gm(2) > 1.5*k2d ) gm(2)=gm(2)-3.0*k2d
          DO i1=0,3*k1d-1
            gm(1)=real(i1)
            IF ( gm(1) > 1.5*k1d ) gm(1)=gm(1)-3.0*k1d
!
!-> use inversion <-> c.c.
!
            ic1 = NINT(gm(1)) ; ic2 = NINT(gm(2)) ; ic3 = NINT(gm(3))
            IF ( gm(3) < 0.0 ) THEN  ! retreive from table icm()
              icc = icm(-ic1,-ic2,-ic3)
              IF (icc.EQ.0) THEN
                write(*,*) ic1,ic2,ic3,icc
                STOP ' error in stepf! '
              ENDIF
              ufft(ic) = ufft(icc)
#ifndef CPP_INVERSION
              bfft(ic) = - bfft(icc)
#endif
              ic=ic+1
              CYCLE 
            ELSE                         ! store number in table icm()
              icm(ic1,ic2,ic3) = ic
              IF (ic1 == im1) icm(-ic1,ic2,ic3) = ic
              IF (ic2 == im2) icm(ic1,-ic2,ic3) = ic
              IF ((ic1 == im1).AND.(ic2 == im2)) icm(-ic1,-ic2,ic3) = ic
            ENDIF
            IF (ic.eq.0) THEN

              ufft(0)=0.0
              bfft(0)=0.0
              DO n=1,ntype
                ufft(0)=ufft(0)+neq(n)*volmts(n)
              ENDDO
              ufft(0)=1.0-ufft(0)*inv_omtil

            ELSE

              CALL dset(3,0.0,g,1)
              CALL CPP_BLAS_sgemv('T',3,3,1.0,bmat,3,gm,1,0.0,g,1)
              g_sqr = CPP_BLAS_sdot(3,g,1,g,1)
              g_abs = sqrt(g_sqr)
              help = fp_omtil/g_sqr
              na=1
#ifdef CPP_INVERSION
              r_c = 0.0
              DO n=1,ntype
                r_phs = 0.0
                DO nn=1,neq(n)
                  th=-tpi*CPP_BLAS_sdot(3,gm,1,taual(1,na),1)
                  r_phs = r_phs + cos(th)
                  na=na+1
                ENDDO
                g_rmt = g_abs * rmt(n)
                r_c=r_c+rmt(n)*(sin(g_rmt)/g_rmt-cos(g_rmt))*r_phs
              ENDDO
              ufft(ic) = help * r_c
#else
              c_c=cmplx(0.0,0.0)
              DO n=1,ntype
                c_phs = cmplx(0.0,0.0)
                DO nn=1,neq(n)
                  th=-tpi*CPP_BLAS_sdot(3,gm,1,taual(1,na),1)
                  c_phs = c_phs + cmplx(cos(th),sin(th))
                  na=na+1
                ENDDO
                g_rmt = g_abs * rmt(n)
                c_c=c_c+rmt(n)*(sin(g_rmt)/g_rmt-cos(g_rmt))*c_phs
              ENDDO
              ufft(ic) = help * real(c_c)
              bfft(ic) = help * aimag(c_c)
#endif
            ENDIF

            IF (((i3.eq.3*k3d/2).or.
     +           (i2.eq.3*k2d/2)).or.
     +           (i1.eq.3*k1d/2)) THEN
                   ufft(ic)=0.0 
                   bfft(ic)=0.0 
            ENDIF
c-odim
            IF (odd%d1) THEN
               IF (ic.LT.9*k1d*k2d .AND. ic.NE.0) THEN
                  gx = (bmat(1,1)*gm(1) + bmat(2,1)*gm(2))
                  gy = (bmat(1,2)*gm(1) + bmat(2,2)*gm(2))
                  gr = sqrt(gx**2 + gy**2)
                  CALL od_cylbes(1,gr*z1,fJ)
                  ufft(ic) = ufft(ic)
     +                 +2*vol*fJ/(gr*z1*omtil)
               END IF
            END IF
c+odim	
            ic=ic+1
          ENDDO
        ENDDO
      ENDDO
c
c --> add film-contributions
c
      IF (film .AND. .NOT.odd%d1) THEN

        ifft2d=9*k1d*k2d
        ufft(0)=ufft(0)+vol*inv_omtil-1.0

        DO i3=1,3*k3d-1
          gm(3)=real(i3)
          IF ( gm(3) > 1.5*k3d ) gm(3)=gm(3)-3.0*k3d
          th=bmat(3,3)*gm(3)*z1
          ufft(i3*ifft2d)=ufft(i3*ifft2d)+vol*inv_omtil*sin(th)/th
        ENDDO

      ELSEIF (odd%d1) THEN
c-odim
         ufft(0) = ufft(0)+vol*inv_omtil-1.0
c+odim

      ENDIF
c
c --> make fft
c
#ifdef CPP_INVERSION
      CALL dset(ifftd,0.0,bfft,1)
#endif
      CALL cfft(ufft,bfft,ifftd,3*k1d,3*k1d,+1)
      CALL cfft(ufft,bfft,ifftd,3*k2d,9*k1d*k2d,+1)
      CALL cfft(ufft,bfft,ifftd,3*k3d,ifftd,+1)

      DEALLOCATE ( bfft , icm )

c--->    store on unit14
      REWIND 14
      WRITE (14) nq3,ifftd
      WRITE (14) (ustep(i),i=1,nq3)
      WRITE (14) (ufft(i),i=0,ifftd-1)

      CLOSE (14)

      END SUBROUTINE stepf
      END MODULE m_stepf
