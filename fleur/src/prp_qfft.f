      MODULE m_prpqfft
      CONTAINS
      SUBROUTINE prp_qfft(
     >                    n2d,n3d,kq1d,kq2d,kq3d,
     >                    ng2,ng3,nstr2,nstr,gmax,bmat,sk2,sk3,kv3,l_ss,
     =                    rkmax,
     <                    mq1d,mq2d,mq3d,nq2_fft,nq3_fft,
     <                    kmxq2_fft,kmxq_fft)
c*********************************************************************
c     This subroutine prepares the necessary variables and checks
c     to calculate the plane wave chargedensity in the interstitial
c     (subroutine pwden.f) by fast fourier transform (FFT).
c     In certain cases rkmax is slightly reajusted to perform
c     a quick FFT based on powers   (2**P) * (3**Q) * (5**R) only.
c
c     dimensions kd(i)d for charge denisity box are checked.
c        s. bluegel, JRCAT, Feb. 97
c
c     subroutine boxdim added to subroutine to treat non-orthogonal
c     lattice systems.
c        s.bluegel, IFF, 18.Nov.97
c*********************************************************************
c
      USE m_ifft, ONLY : ifft235

      IMPLICIT NONE
c
c---> Arguments
c
      INTEGER, INTENT (IN)  :: n2d,n3d,kq1d,kq2d,kq3d,ng2,ng3
      LOGICAL, INTENT (IN)  :: l_ss
      REAL,    INTENT (IN)  :: gmax
      REAL,    INTENT (INOUT) :: rkmax
      INTEGER, INTENT (OUT) :: mq1d,mq2d,mq3d
      INTEGER, INTENT (OUT) :: nq2_fft,nq3_fft,kmxq2_fft,kmxq_fft
!
      INTEGER, INTENT (IN) :: kv3(3,n3d),nstr2(n2d),nstr(n3d)
      REAL,    INTENT (IN) :: bmat(3,3),sk2(n2d),sk3(n3d)
c
c---> local variables
c
      INTEGER ksfft,mq1,mq2,mq3,istr,iofile
      REAL    arltv1,arltv2,arltv3,rknew
c
c---> external subroutines ..
c
      EXTERNAL boxdim
c
c---> intrinsic functions
c    
      INTRINSIC int,min
c
c---> data statement
c
      REAL gmaxp
      DATA gmaxp /2.0/
c
C------->          BACKGROUND       
c
c        Determine the limits  for the charge density fft-box
c        The charge density is the  square of the wavefunction.
c        Since the largest G-vector of the wavefunction is
c        given by |G + k| < rkmax, the largest G -vector gmaxq
c        contributing the charge-density is exactly :
c        gmaxq = gmaxp*rkmax,  with gmaxp =2. 
c        Therefore the box dimension must be  determined by :
c        m(i) >= gmaxp * rkmax*a(i)/(2*pi), where a(i) is the 
c        magnitude of the ith basis vector. (this is only true
c        for <a(i)|a(j)> = 0 for i .ne. j., otherwise see boxdim)
c        REMEMBER: gmax used through out the FLAPW code can be
c        larger than gmaxp*rkmax. 
c
C
C------->          ABBREVIATIONS
C
C   ksfft=(0,1) : KEY OF SELECTING FFT-PRDOGRAM AND RADIX-TYPE
C                      0  PROGRAM, RADIX-2 ONLY
C                      1  PROGRAM, RADIX-2, RADIX-3,RADIX-5
c   iofile      : device number for in and output
c   gmax        : actually used gmax troughtout the FLAPW code
c   gmaxq       : cut-off wavevector for charge density
c   rkmax       : cut-off for |g+k|
c   gmaxp       : gmaxp = gmaxq/rkmax, ideal: gmaxp=2
c   kq1,2,3d    : dimensions of the charge density fft-box
c                 ( in positive fft domain )
c   mq1,2,3     : actual size of the charge density fft-box
c                 ( in positive fft domain )
c   nq3_fft     : number of stars in the charge density sphere
c   kmxq_fft    : number of g-vectors in the charge density sphere
c   arltv(i)    : length of reciprical lattice vector along
c                 direction (i)
C
C---> Determine dimensions of fft-box of size mq1, mq2, mq3,
c     for which |G(mq1,mq2,mq3)| < Gmaxp*Kmax
C
      CALL boxdim(
     >            bmat,
     <            arltv1,arltv2,arltv3)
c
      mq1 = int( gmaxp*rkmax/arltv1 ) + 1
      mq2 = int( gmaxp*rkmax/arltv2 ) + 1
      mq3 = int( gmaxp*rkmax/arltv3 ) + 1

c---> add + 1 in spin spiral calculation, to make sure that all G's are 
c---> still within the FFT-box after being shifted by the spin spiral
c---> q-vector.
      IF (l_ss) THEN
         mq1 = mq1 + 1
         mq2 = mq2 + 1
         mq3 = mq3 + 1
      ENDIF         
c
c----> fft done in positive domain 
c      (Remember number of grid points not 2*m+1 because of cyclic 
c       boundary condition)
c
      mq1 = mq1 + mq1
      mq2 = mq2 + mq2
      mq3 = mq3 + mq3
c
c---> fft's are usually fastest for low primes
c     (restrict kqid to: kwid=  (2**P) * (3**Q) * (5**R)
c
      iofile = 6
      ksfft  = 1
      mq1d = ifft235 (iofile,ksfft,mq1,gmaxp)
      mq2d = ifft235 (iofile,ksfft,mq2,gmaxp)
      mq3d = ifft235 (iofile,ksfft,mq3,gmaxp)
c+gb
!      mq1d = mq1
!      mq2d = mq2
!      mq3d = mq3
c
c---> For mq = 2**P, FFT is very fast. If mq very close to 2**P
c     choose this, even 2**P < mq . Therefore:
c     factor 0.5 added by S.B. 21.Aug.97 
c
      rknew = min( real( mq1d )*arltv1, real( mq2d )*arltv2, 
     +             real( mq3d )*arltv3)
      rknew = 0.5 * rknew / gmaxp
      IF (rknew.LT.rkmax) THEN
         WRITE (6,'('' rkmax and true gmax recalculated '')')
         WRITE (6,2100) rkmax, rknew, rknew*rknew
         WRITE (6,2200) gmaxp*rknew, gmaxp*rknew*gmaxp*rknew
         WRITE (16,'('' rkmax and true gmax recalculated '')')
         WRITE (16,2100) rkmax, rknew, rknew*rknew
         WRITE (16,2200) gmaxp*rknew, gmaxp*rknew*gmaxp*rknew
         rkmax = rknew
      ENDIF
c
c-----> gmax => gmaxp*rkmax 
c       (otherwise too few elements in arrays defined in strng1)
c
      IF (rkmax+rkmax.GT.gmax) THEN
         WRITE (6,'('' gmax must be at least 2*rkmax'')')
         WRITE (6,'('' increase gmax , or reduce rkmax'')')
         WRITE (6,'('' rkmax ='',f10.3,''  gmax ='',f10.3)') rkmax,gmax
         STOP 'rkmax,gmax'
      ENDIF
c
c     mq1 = mq1d
c     mq2 = mq2d
c     mq3 = mq3d
c
c------> check dimensions of FFT chargedensity box used in pwden.f
c
       IF ( mq1d.gt.kq1d .OR. mq2d.gt.kq2d .OR. mq3d.gt.kq3d) THEN
          WRITE ( 6,'('' box dim. for FFT too small'')')
          WRITE ( 6,'('' mq1d,mq2d,mq3d,kq1d,kq2d,kq3d '',6i5)')
     +                   mq1d,mq2d,mq3d,kq1d,kq2d,kq3d
          WRITE (16,'('' box dim. for FFT too small'')')
          WRITE (16,'('' mq1d,mq2d,mq3d,kq1d,kq2d,kq3d '',6i5)')
     +                 mq1d,mq2d,mq3d,kq1d,kq2d,kq3d
          STOP 'prp_qfft'
       ENDIF
c
c-----> how many stars are in charge density sphere?
c       assume stars are ordered according to length
c
c---> 3d stars
      nq3_fft = 0
      DO istr = 1 , ng3
         IF ( sk3(istr).LE.gmaxp*rkmax ) THEN
            nq3_fft = istr
         ENDIF
      ENDDO
c---> 2d stars
      nq2_fft = 0
      DO istr = 1 , ng2
         IF ( sk2(istr).LE.gmaxp*rkmax ) THEN
            nq2_fft = istr
         ENDIF
      ENDDO
c
      IF ( nq3_fft.EQ.0 ) THEN
         WRITE(6,'('' presumably ng3 too small '')')
         WRITE(6,'('' sk3max, gmaxp*rkmax '', 2f10.6)') 
     +                sk3(ng3),gmaxp*rkmax
         WRITE(16,'('' presumably ng3 too small '')')
         WRITE(16,'('' sk3max, gmaxp*rkmax '', 2f10.6)')
     +                sk3(ng3),gmaxp*rkmax
         STOP 'prp_qfft'
      ENDIF
c
      IF ( nq3_fft.GT.n3d ) THEN
         WRITE(6,'('' nq3_fft > n3d '')')
         WRITE(6,'('' nq3_fft, n3d '',2i10)') nq3_fft, n3d
         WRITE(16,'('' nq3_fft > n3d '')')
         WRITE(16,'('' nq3_fft, n3d '',2i10)') nq3_fft, n3d
         STOP 'prp_qfft'
      ENDIF
c
c-----> check that all nq3_fft stars fit into the charge density FFT-box
c
      DO istr = 1 , nq3_fft
        IF ( ( 2*kv3(1,istr).GT.mq1d ) .OR.
     +       ( 2*kv3(2,istr).GT.mq2d ) .OR.
     +       ( 2*kv3(3,istr).GT.mq3d ) ) THEN
          WRITE (6,'('' not all nq3_fft stars in chg. den. FFT box'')')
          WRITE (6,'('' inconsistency in def.s see also strgn1'')')
          WRITE (6,'('' mq1d,mq2d,mq3d,kv1,kv2,kv3 '',6i5)')
     +                  mq1d,mq2d,mq3d,2*kv3(1,istr),2*kv3(2,istr),
     +                                 2*kv3(3,istr)
          WRITE (16,'('' not all nq3_fft stars in chg. den. FFT box'')')
          WRITE (16,'('' inconsistency in def.s see also strgn1'')')
          WRITE (16,'('' mq1d,mq2d,mq3d,kv1,kv2,kv3 '',6i5)')
     +                  mq1d,mq2d,mq3d,2*kv3(1,istr),2*kv3(2,istr),
     +                                 2*kv3(3,istr)
          STOP 'prp_qfft'
        ENDIF
      ENDDO
c
c-----> how many G-vectors belong to these nq3_fft stars                    
c
c---> 3d vectors
      kmxq_fft = 0
      DO istr = 1 , nq3_fft
         kmxq_fft = kmxq_fft + nstr(istr)
      ENDDO
      IF ( kmxq_fft .GT. kq1d*kq2d*kq3d ) THEN
         WRITE (6,'('' array dimensions in later subroutines too'',
     +             '' small'',2i10)') kmxq_fft,kq1d*kq2d*kq3d
      ENDIF
c---> 2d vectors
      kmxq2_fft = 0
      DO istr = 1 , nq2_fft
         kmxq2_fft = kmxq2_fft + nstr2(istr)
      ENDDO
      IF ( kmxq2_fft .GT. kq1d*kq2d ) THEN
         WRITE (6,'('' array dimensions in later subroutines too'',
     +             '' small'',2i10)') kmxq2_fft,kq1d*kq2d
      ENDIF

c
 2100 FORMAT (/,1x,'old rkmax   =',f10.5, '(a.u.)**(-1) ==>  new rkmax '
     >      ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new E_cut(wf)   =',
     >            f10.5,' Ry')
 2200 FORMAT (/,1x,'true gmax   =',f10.5, '(a.u.)**(-1)',
     >       /,t38,'==>  new E_cut(chg)  =', f10.5,' Ry')

      END SUBROUTINE prp_qfft
      END MODULE m_prpqfft
