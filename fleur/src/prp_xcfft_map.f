      MODULE m_prpxcfftmap
c*********************************************************************
c     this subroutine prepares the pointer which identifies a
c     threedimensional g-vector in the positive domain of the
c     xc (=charge density) fft-box in order to map a 3-d g-vector
c     onto stars in case of the backtransform for fft of the
c     charge density. correspondes  to igfft(*,2)
c     Further it sets up the x,y, and z component of the 3-dimensional
c     g-vector in the original domain of all g-vectors used for fft.
c     it is independent of spin and k-point.
c     pointer is built up when ever the chargedensity is calculated
c     in order to save memory
c
c        s. bluegel, IFF, Aug. 97   
c*********************************************************************
      CONTAINS
      SUBROUTINE prp_xcfft_map(
     >                         n3d,kxc1d,kxc2d,kxc3d,nop,
     >                         kv3,nop2,mrot,bmat,
     >      kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft,
     <                         igxc_fft,gxc_fft)
c
      IMPLICIT NONE
c
      INTEGER, INTENT (IN) :: n3d,kxc1d,kxc2d,kxc3d,nop
c
      INTEGER nxc3_fft,kmxxc_fft,kv3(3,n3d),nop2,mrot(3,3,nop)
      INTEGER kxc1_fft,kxc2_fft,kxc3_fft,igxc_fft(0:kxc1d*kxc2d*kxc3d-1)
      REAL    bmat(3,3),gxc_fft(0:kxc1d*kxc2d*kxc3d-1,3)
c
c---> local variables
c
      LOGICAL new
      INTEGER istr,iop,iopm1,il,im,in,kidx,iv1d,ifftq1,ifftq2
      INTEGER nop_local,norm,kr(3,nop)

c------->          abbreviations
c
c     kxc1d  : dimension of the charge density fft box in the pos. domain
c     kxc2d  : defined in dimens.f program (subroutine apws). 1,2,3 indic
c     kxc3d  ; a_1, a_2, a_3 directions.
c     kq(i) : i=1,2,3 actual length of the fft-box for which fft is done
c     nstr  : number of members (arms) of reciprocal lattice (g) vector
c             of each star.
c     nxc3_fft: number of stars in the  charge density  fft-box
c     kmxxc_fft: number of g-vectors forming the nxc3_fft stars in the
c               charge density sphere
c     gxc_fft : contains x,y,z components of g-vectors contributing to FFT.
c
c-----> prepare pointer which identifies a threedimensional g-vector
c       in the positive domain of the charge density fft-box.
c       correspondes  to igfft(*,2)
c
      kidx    = 0
      ifftq1  = kxc1_fft
      ifftq2  = kxc1_fft*kxc2_fft
c
      DO istr = 1 , nxc3_fft
c
         nop_local=nop
         IF (kv3(3,istr).EQ.0) nop_local=nop2
c
         DO iop = 1,nop_local
            kr(1,iop) = kv3(1,istr)*mrot(1,1,iop)
     +                + kv3(2,istr)*mrot(2,1,iop)
     +                + kv3(3,istr)*mrot(3,1,iop)
            kr(2,iop) = kv3(1,istr)*mrot(1,2,iop)
     +                + kv3(2,istr)*mrot(2,2,iop)
     +                + kv3(3,istr)*mrot(3,2,iop)
            kr(3,iop) = kv3(1,istr)*mrot(1,3,iop)
     +                + kv3(2,istr)*mrot(2,3,iop)
     +                + kv3(3,istr)*mrot(3,3,iop)
         ENDDO
c
         DO iop = 1 , nop_local
            new=.true.
            DO iopm1 = 1 , iop - 1
              norm=(kr(1,iop)-kr(1,iopm1))**2 +
     +             (kr(2,iop)-kr(2,iopm1))**2 +
     +             (kr(3,iop)-kr(3,iopm1))**2
              if(norm.eq.0) new=.false.
            ENDDO

            IF (new) THEN
              il=kr(1,iop)
              im=kr(2,iop)
              in=kr(3,iop)
              gxc_fft(kidx,1) = bmat(1,1)*il+bmat(2,1)*im+bmat(3,1)*in 
              gxc_fft(kidx,2) = bmat(1,2)*il+bmat(2,2)*im+bmat(3,2)*in
              gxc_fft(kidx,3) = bmat(1,3)*il+bmat(2,3)*im+bmat(3,3)*in
              IF (il.LT.0) il=il+kxc1_fft
              IF (im.LT.0) im=im+kxc2_fft
              IF (in.LT.0) in=in+kxc3_fft
              iv1d = in*ifftq2 + im*ifftq1 + il
              igxc_fft(kidx)=iv1d
              kidx=kidx+1
            ENDIF
         ENDDO
      ENDDO
c
      IF (kidx .NE. kmxxc_fft) THEN
         WRITE (6,'('' something wrong with kmxxc_fft or nxc3_fft'')')
         WRITE (6,'('' kmxxc_fft, acutal kidx '',2i5)')
     +                kmxxc_fft, kidx
         STOP 'prp_xcfft_map: kidx .ne. kmxxc_fft'
      ENDIF
c
      END SUBROUTINE prp_xcfft_map
      END MODULE m_prpxcfftmap
