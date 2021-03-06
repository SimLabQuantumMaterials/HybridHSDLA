      MODULE m_fft3dxc
      CONTAINS
      SUBROUTINE fft3dxc(
     x                 afft,bfft,fg3,
     >                 k1d,k2d,k3d,ng3,kimax,isn,
     >                 igfft1,igfft2,pgfft,nstr)

************************************************************
*                                                          *
* interface for fg3(star) -- fft --> (a,b)fft (r) (isn=+1) *
*         or (a,b)fft (r) -- fft --> fg3(star)    (isn=-1) *
*                                                          *
* dimension of (a,b)fft is (k1d x k2d x k3d)               *
* afft and bfft contain the real/imaginary part of the fft *
* igfft1(i)  is the pointer from the g-sphere to stars     *
* igfft2(i)  is the pointer from the g-sphere to fft-grid  *
* pgfft(i)   contains the phases of the g-vectors of sph.  *
*                                                          *
************************************************************

      USE m_set, ONLY : dset
      IMPLICIT NONE

      INTEGER k1d,k2d,k3d,ng3,kimax,isn
      INTEGER igfft1(0:k1d*k2d*k3d-1),igfft2(0:k1d*k2d*k3d-1)
      INTEGER nstr(ng3)
      REAL    pgfft(0:k1d*k2d*k3d-1)
      REAL    afft(0:k1d*k2d*k3d-1),bfft(0:k1d*k2d*k3d-1)
      COMPLEX fg3(ng3)

      INTEGER i,ifftd
      REAL scale,zero

      ifftd=k1d*k2d*k3d
      zero=0.0

      IF( isn.GT.0) THEN
c
c  ---> put stars onto the fft-grid
c
        CALL dset(ifftd,zero,afft,1)
        CALL dset(ifftd,zero,bfft,1)
        DO i=0,kimax-1
          afft(igfft2(i))=real(fg3(igfft1(i)))  * pgfft(i)
          bfft(igfft2(i))=aimag(fg3(igfft1(i))) * pgfft(i)
        ENDDO
      ENDIF

c---> now do the fft (isn=+1 : g -> r ; isn=-1 : r -> g)

      CALL cfft(afft,bfft,ifftd,k1d,k1d,isn)
      CALL cfft(afft,bfft,ifftd,k2d,k1d*k2d,isn)
      CALL cfft(afft,bfft,ifftd,k3d,ifftd,isn)

      IF (isn.LT.0) THEN
c
c  ---> collect stars from the fft-grid
c
        DO i=1,ng3
          fg3(i) = cmplx(0.0,0.0)
        ENDDO
        scale=1.0/ifftd
        DO i=0,kimax-1
          fg3(igfft1(i))=fg3(igfft1(i))+pgfft(i)*
     +                cmplx(afft(igfft2(i)),bfft(igfft2(i)))
        ENDDO
        DO i=1,ng3
          fg3(i)=scale*fg3(i)/nstr(i)
        ENDDO
      ENDIF

      END SUBROUTINE fft3dxc
      END MODULE m_fft3dxc
