      MODULE m_fft2d
      CONTAINS
      SUBROUTINE fft2d(
     >                 k1d,k2d,n2d,
     X                 afft2,bfft2,
     X                 fg,fgi,fgxy,
     >                 stride,ng2,kimax2,isn,
     >                 igfft2,pgfft2,nstr2)

*************************************************************
*                                                           *
* interface for fg2(star) -- FFT --> gfft (r)     (isn=+1)  *
*            or gfft (r)  -- FFT --> fg2(star)    (isn=-1)  *
*                                                           *
* dimension of gfft2 is (3*k1d x 3*k2d)                     *
* afft/bfft contain the real/imaginary part of gfft         *
* igfft2(i,1) is the pointer from the G-sphere to stars     *
* igfft2(i,2) is the pointer from the G-sphere to fft-grid  *
* pgfft2(i)   contains the phases of the G-vectors of sph.  *
*                                                           *
*************************************************************

#include"cpp_double.h"
      USE m_set, ONLY : dset
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: k1d,k2d,n2d
      INTEGER, INTENT (IN) :: ng2,kimax2,isn,stride
      REAL                 :: fg,fgi

      INTEGER,INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      INTEGER,INTENT (IN) :: nstr2(ng2)
      REAL,   INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL,   INTENT (INOUT):: afft2(0:9*k1d*k2d-1),bfft2(0:9*k1d*k2d-1)
      COMPLEX               :: fgxy(stride,ng2-1)

C... local variables
      
      INTEGER i,ifftd2
      REAL scale
      COMPLEX fg2(n2d)

      ifftd2=9*k1d*k2d
c
      IF (isn.gt.0) THEN
c
C  ---> put stars onto the fft-grid 
c
        fg2(1) = cmplx(fg,fgi)
        CALL CPP_BLAS_ccopy(ng2-1,fgxy,stride,fg2(2),1)
        CALL dset(ifftd2,0.0,afft2,1)
        CALL dset(ifftd2,0.0,bfft2,1)
        DO i=0,kimax2-1
          afft2(igfft2(i,2))=real(fg2(igfft2(i,1)))*
     +                                       pgfft2(i)
          bfft2(igfft2(i,2))=aimag(fg2(igfft2(i,1)))*
     +                                       pgfft2(i)
        ENDDO
      ENDIF

C---> now do the fft (isn=+1 : G -> r ; isn=-1 : r -> G)

      CALL cfft(afft2,bfft2,ifftd2,3*k1d,3*k1d,isn)
      CALL cfft(afft2,bfft2,ifftd2,3*k2d,ifftd2,isn)

      IF (isn.lt.0) THEN
c
C  ---> collect stars from the fft-grid
c
        DO i=1,ng2
          fg2(i) = cmplx(0.0,0.0)
        ENDDO
        scale=1.0/ifftd2
        DO i=0,kimax2-1
          fg2(igfft2(i,1))=fg2(igfft2(i,1))+pgfft2(i)*
     +                 cmplx(afft2(igfft2(i,2)),bfft2(igfft2(i,2)))
        ENDDO
        fg=scale*real(fg2(1))/nstr2(1)
        fgi=scale*aimag(fg2(1))/nstr2(1)
        DO i=2,ng2
          fgxy(1,i-1)=scale*fg2(i)/nstr2(i)
        ENDDO
      ENDIF

      END SUBROUTINE fft2d
      END MODULE m_fft2d
