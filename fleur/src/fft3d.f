      MODULE m_fft3d
      CONTAINS
      SUBROUTINE fft3d(
     X                 afft,bfft,fg3,
     >                 k1d,k2d,k3d,ng3,kimax,isn,
     >                 igfft,pgfft,nstr)

************************************************************
*                                                          *
* interface for fg3(star) -- FFT --> (a,b)fft (r) (isn=+1) *
*         or (a,b)fft (r) -- FFT --> fg3(star)    (isn=-1) *
*                                                          *
* dimension of (a,b)fft is (3*k1d x 3*k2d x 3*k3d)         *
* afft and bfft contain the real/imaginary part of the FFT *
* igfft(i,1) is the pointer from the G-sphere to stars     *
* igfft(i,2) is the pointer from the G-sphere to fft-grid  *
* pgfft(i)   contains the phases of the G-vectors of sph.  *
*                                                          *
************************************************************

      USE m_set, ONLY : dset
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: k1d,k2d,k3d,ng3,kimax,isn
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      INTEGER, INTENT (IN) :: nstr(ng3)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (INOUT) :: afft(0:27*k1d*k2d*k3d-1)
      REAL,    INTENT (INOUT) :: bfft(0:27*k1d*k2d*k3d-1)
      COMPLEX                 :: fg3(ng3)

      INTEGER i,ifftd
      REAL scale,zero

      ifftd=27*k1d*k2d*k3d
      zero=0.0

      IF (isn.GT.0) THEN
c
C  ---> put stars onto the fft-grid 
c
        CALL dset(ifftd,zero,afft,1)
        CALL dset(ifftd,zero,bfft,1)
        DO i=0,kimax-1
          afft(igfft(i,2))=real(fg3(igfft(i,1)))*
     +                                   pgfft(i)
          bfft(igfft(i,2))=aimag(fg3(igfft(i,1)))*
     +                                   pgfft(i)
        ENDDO
      ENDIF

C---> now do the fft (isn=+1 : G -> r ; isn=-1 : r -> G)

      CALL cfft(afft,bfft,ifftd,3*k1d,3*k1d,isn)
      CALL cfft(afft,bfft,ifftd,3*k2d,9*k1d*k2d,isn)
      CALL cfft(afft,bfft,ifftd,3*k3d,ifftd,isn)

      IF (isn.LT.0) THEN
c
C  ---> collect stars from the fft-grid
c
        DO i=1,ng3
          fg3(i) = cmplx(0.0,0.0)
        ENDDO
        DO i=0,kimax-1
          fg3(igfft(i,1))=fg3(igfft(i,1))+  pgfft(i)*
     +                cmplx(afft(igfft(i,2)),bfft(igfft(i,2)))
        ENDDO
        scale=1.0/ifftd
        DO i=1,ng3
          fg3(i)=scale*fg3(i)/nstr(i)
        ENDDO
      ENDIF

      END SUBROUTINE fft3d
      END MODULE m_fft3d
