      MODULE m_prpqfftmap
      CONTAINS
      SUBROUTINE prp_qfft_map(
     >                        n3d,kq1d,kq2d,kq3d,nop,
     >                        kv3,nop2,mrot,film,
     >                        kq1,kq2,kq3,
     >                        nq3_fft,kmxq2_fft,kmxq_fft,
     <                        igq2_fft,igq_fft)
c*********************************************************************
c     This subroutine prepares the pointer which identifies a 
c     threedimensional g-vector in the positive domain of the 
c     charge density fft-box in order to map a 3-d g-vector
c     onto stars in case of the backtransform for fft of the 
c     charge density. correspondes  to igfft(*,2)     
c     it is independent of spin and k-point. 
c     pointer is built up when ever the chargedensity is calculated
c     in order to save memory
c
c        s. bluegel, JRCAT, Feb. 97
c*********************************************************************
c
      IMPLICIT NONE
c
      INTEGER, INTENT (IN) :: n3d,kq1d,kq2d,kq3d,nop
c
      INTEGER nq3_fft,kmxq2_fft,kmxq_fft,kv3(3,n3d)
      INTEGER nop2,mrot(3,3,nop)
      INTEGER kq1,kq2,kq3
      INTEGER igq2_fft(0:kq1d*kq2d-1),igq_fft(0:kq1d*kq2d*kq3d-1)
      LOGICAL film
c
c---> local variables
c
      LOGICAL new
      INTEGER istr,iop,iopm1,il,im,in,kid2x,kidx,iv1d,ifftq1,ifftq2
      INTEGER norm,kr(3,nop),nop_local

C------->          ABBREVIATIONS
C
c     kq1d  : dimension of the charge density FFT box in the pos. domain
c     kq2d  : defined in dimens.f program (subroutine apws). 1,2,3 indicate
c     kq3d  ; a_1, a_2, a_3 directions.
c     kq(i) : i=1,2,3 actual length of the fft-box for which FFT is done.
c     nstr  : number of members (arms) of reciprocal lattice (g) vector
c             of each star.
c     nq3_fft: number of stars in the  charge density  FFT-box
c     kmxq_fft: number of g-vectors forming the nq3_fft stars in the
c               charge density sphere
c
c-----> prepare pointer which identifies a threedimensional g-vector
c       in the positive domain of the charge density fft-box.
c       correspondes  to igfft(*,2)     
c
      kidx    = 0
      kid2x   = 0
      ifftq1  = kq1
      ifftq2  = kq1*kq2
c
      DO istr = 1 , nq3_fft
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
     +                (kr(2,iop)-kr(2,iopm1))**2 +
     +                (kr(3,iop)-kr(3,iopm1))**2
                 IF (norm.EQ.0) new=.false.
               ENDDO
 
               IF (new) THEN
                 il=kr(1,iop)
                 im=kr(2,iop)
                 in=kr(3,iop)
                 if(il.lt.0) il=il+kq1
                 if(im.lt.0) im=im+kq2
                 if(in.lt.0) in=in+kq3
                 iv1d = in*ifftq2 + im*ifftq1 + il
                 igq_fft(kidx)=iv1d 
                 kidx=kidx+1
                 IF (film.AND.(kv3(3,istr).EQ.0)) THEN
                    iv1d = im*ifftq1 + il
                    igq2_fft(kid2x)=iv1d 
                    kid2x=kid2x+1
                 ENDIF
               ENDIF
         ENDDO
c
      ENDDO
c
      IF (kidx .NE. kmxq_fft) THEN
         WRITE (6,'('' something wrong with kmxq_fft or nq3_fft'')')
         WRITE (6,'('' kmxq_fft, acutal kidx '',2i5)') 
     +                kmxq_fft, kidx
         STOP 'prp_qfft_map'
      ENDIF
      IF (film.AND.(kid2x .NE. kmxq2_fft).AND.(kmxq2_fft.NE.0)) THEN
         WRITE (6,'('' something wrong with kmxq2_fft or nq2_fft'')')
         WRITE (6,'('' kmxq2_fft, acutal kid2x '',2i5)') 
     +                kmxq2_fft, kid2x
         STOP 'prp_qfft_map'
      ENDIF         

      END SUBROUTINE prp_qfft_map
      END MODULE m_prpqfftmap
