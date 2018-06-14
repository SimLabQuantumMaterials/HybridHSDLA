      MODULE m_convol
      CONTAINS
      SUBROUTINE convol(
     >                  k1d,k2d,k3d,n3d, 
     <                  fg3,
     >                  ag3,ng3,
     =                  kimax,igfft,pgfft,ufft)

************************************************************
*                                                          *
* calculate f(G) = \sum_G' U(G' - G) a(G')                 *
*                                                          *
*       ag3(star) -- FFT --> gfft(r,1)                     *
*                            gfft(r,1)=gfft(r,1) * U (r)   *
*       fg3(star) <- FFT --- gfft(r,1)                     *
*                                                          *
* dimension of gfft is (3*k1d x 3*k2d x 3*k3d)             *
*                                                          *
************************************************************

      USE m_fft3d
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d
      INTEGER, INTENT (IN) :: ng3,kimax

      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
      COMPLEX, INTENT (IN) :: ag3(n3d)
      COMPLEX, INTENT (OUT):: fg3(n3d)

      INTEGER i,ifftd
      INTEGER nstr(n3d)
      REAL, ALLOCATABLE :: gfft(:,:)

      ifftd=27*k1d*k2d*k3d
      DO i=1,ng3
        nstr(i)=1
      ENDDO
      ALLOCATE (gfft(0:ifftd-1,2))

      CALL fft3d(
     <           gfft(0,1),gfft(0,2),
     >           ag3,
     >           k1d,k2d,k3d,ng3,kimax,+1,
     >           igfft,pgfft,nstr) 

      DO i=0,ifftd-1
        gfft(i,:)=gfft(i,:)*ufft(i)
      ENDDO

      CALL fft3d(
     >           gfft(0,1),gfft(0,2),
     <           fg3,
     >           k1d,k2d,k3d,ng3,kimax,-1,
     >           igfft,pgfft,nstr) 

      DEALLOCATE (gfft)

      END SUBROUTINE convol
      END MODULE m_convol
