      MODULE m_rfft
      CONTAINS
      SUBROUTINE rfft(
     >                isn,n1d,n2d,n3d,n1,n2,n3,
     >                nw1,nw2,nw3,wsave,b,
     X                a)
c **********************************************************************
c
c isn = +1 : FFT a real, 3dim array "a" of dimensions n1*n2*n3 to a complex
c            array of size n1*n2*(n3+1)/2 [if n3 is odd] or n1*n2*(n3/2+1)
c            [n3 is even]. 
c isn = -1 : back-FFT of a complex array "a" of size n1*n2*(n3+1)/2 [odd] 
c            or n1*n2*(n3/2+1) [even] to a real, 3dim array 
c
c the actual array is assumed to be located between 
c                                  1 ... nw1+1    &     n1-nw1+1 ... n1
c                                  1 ... nw2+1    &     n2-nw2+1 ... n2
c                                  1 ... nw3+1    &     n3-nw3+1 ... n3
c and padded with zeros in between.
c if nw1 >= (n1-1)/2, no padding is assumed (-> aliasing errors!)
c
c                                                  G.Bihlmayer (UniWien)
c         
c **********************************************************************
#include"cpp_double.h"
      USE m_set, ONLY : dset
      IMPLICIT NONE

      INTEGER n1d,n2d,n3d,n1,n2,n3,nw1,nw2,nw3,isn
      REAL a(n1d,n2d,0:n3d),b(n1d,n2d,n3d),wsave(n3d+15)

      INTEGER i1,i2,i3,nup
      REAL factor
      LOGICAL l_nopad
c
c a ... array for FFT
c b ... work array
c wsave ... stores tables for r-2-c FFT
c n1,n2,n3 ... dimensions to be transformed
c nw1,nw2,nw3 ... actual dimensions of a before FFT
c n1d,n2d,n3d ... dimensions of a,b 
c
      EXTERNAL vrfftf,vrfftb,vrffti,cfft,CPP_BLAS_sscal
c
c check for input errors
c
      IF ((isn.ne.-1).AND.(isn.ne.1)) STOP 'dfft: choose isn = +/- 1'
      IF ((n1d.lt.n1).OR.(n2d.lt.n2).OR.(n3d.lt.n3)) THEN
        WRITE (6,*) 'n1d,n2d,n3d =',n1d,n2d,n3d
        WRITE (6,*) 'n1 ,n2 ,n3  =',n1 ,n2 ,n3 
        STOP 'dfft: n(i) > n(i)d'
      ENDIF
      IF ((n1.le.2*nw1+1).OR.
     +    (n2.le.2*nw2+1).OR.
     +    (n3.le.2*nw3+1)) THEN
c        WRITE (6,*) 'n1 ,n2 ,n3  =',n1 ,n2 ,n3 
c        WRITE (6,*) 'nw1,nw2,nw3 =',nw1,nw2,nw3
        l_nopad=.true.
      ELSE
        l_nopad=.false.
      ENDIF
c
c ******************** f o r w a r d - t r a n s f o r m *******************     
c      
      IF (isn.eq.1) THEN
c
c array a is assumed to be zero from (1,1,0) to (n1d,n2d,0) and the array 
c to be FFT'd starts at (1,1,1) as n1*n2*n3 real numbers. 
c first transform n1*n2 real sequences of lenghth n3 to n3/2 complex values
c
        CALL vrffti(n3,wsave) 
        IF (l_nopad) THEN
          CALL vrfftf(n1*n2,n3,a(1,1,1),b,n1d*n2d,wsave)
        ELSE
          DO i2=1,nw2+1
            CALL vrfftf(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftf(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
          DO i2=n2-nw2+1,n2
            CALL vrfftf(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftf(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
        ENDIF
c
c now we have the FFT'd array stored as described in vrfftf 
c (mixed real & compex data)
c remove the norm 1/sqrt(n3) (to be compatible with cfft)
c
        factor = sqrt(1.0*n3)
        CALL CPP_BLAS_sscal(n1d*n2d*n3,factor,a(1,1,1),1)
c
c now, the real part of f(0) has to be moved to a(n1,n2,0) to get a purely
c complex array starting at a(1,1,0)
c
        DO i1=1,n1
          DO i2=1,n2
            a(i1,i2,0)=a(i1,i2,1) 
            a(i1,i2,1)=0.0 
          ENDDO
        ENDDO
c
      ENDIF
c
c ******************** g e n e r a l  p a r t *******************
c
c now perform n2*n3/2 and n1*n3/2 complex FFT's; a is assumed to be
c complex and starting at (1,1,0)
c
      IF (ABS((n3/2.)-NINT(n3/2.)).gt.0.1) THEN
        nup = n3
      ELSE
        nup = n3+1
        IF (n3+1.gt.n3d) STOP 'n3 even & n3+1 > n3d'
        CALL dset(n1*n2,0.0,a(1,1,n3+1),1)
      ENDIF

      IF (l_nopad) THEN
        DO i3=1,nup,2
          CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n1,n1,isn)
          CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n2,n1*n2,isn)
        ENDDO
      ELSE
        DO i3=1,nup,2
          CALL cfft(a(1,1,i3-1),a(1,1,i3),(nw2+1)*n1,n1,n1,isn)
          CALL cfft(a(1,n2-nw2+1,i3-1),a(1,n2-nw2+1,i3),
     +                                        nw2*n1,n1,n1,isn)
          CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n2,n1*n2,isn)
        ENDDO
      ENDIF

c
c ******************** b a c k w a r d - t r a n s f o r m *******************     
c      
      IF (isn.eq.-1) THEN
c
c the real part of f(0) has to be moved to a(n1,n2,1) for a correct
c setup for vrfftb (see comments therein)
c
        DO i1=1,n1
          DO i2=1,n2
            a(i1,i2,1)=a(i1,i2,0)
            a(i1,i2,0)=0.0
          ENDDO
        ENDDO
c
c transform n1*n2 mixed real and complex sequences of lenghth n3/2
c to n3 real values
c
        CALL vrffti(n3,wsave)
        IF (l_nopad) THEN
          CALL vrfftb(n1*n2,n3,a(1,1,1),b,n1d*n2d,wsave)
        ELSE
          DO i2=1,nw2+1
            CALL vrfftb(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftb(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
          DO i2=n2-nw2+1,n2
            CALL vrfftb(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftb(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
        ENDIF
c
c remove the norm 1/sqrt(n3) (compatibility with cfft)
c
        factor = sqrt(1.0*n3)
        CALL CPP_BLAS_sscal(n1d*n2d*n3,factor,a(1,1,1),1)
c
      ENDIF

      END SUBROUTINE rfft
      END MODULE m_rfft
