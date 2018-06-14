      MODULE m_prpxcfft
      CONTAINS
      SUBROUTINE prp_xcfft(
     >                     n3d,kxc1d,kxc2d,kxc3d,
     >                     ng3,nstr,gmax,rkmax,bmat,sk3,kv3,
     =                     gmaxxc,
     <                  kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft)
c*********************************************************************
c     this subroutine prepares the necessary variables and checks
c     to calculate the xc-potential and energy in the interstitial
c     (subroutine visxc(g).f) by fast fourier transform (fft).
c     in certain cases gmaxxc is slightly reajusted to perform
c     a quick fft based on powers   (2**p) * (3**q) * (5**r) only.
c
c     dimensions kd(i)d for charge denisity box are checked.
c        s. bluegel, iff , oct. 97
c     subroutine boxdim added to subroutine to treat non-othogonal
c     lattice systems
c        s.bluegel, IFF, 18.Nov.97
c*********************************************************************
c
      USE m_ifft, ONLY : ifft235

      IMPLICIT NONE
c
c
c---> arguments
c
      INTEGER, INTENT (IN) :: n3d,kxc1d,kxc2d,kxc3d
      INTEGER, INTENT (IN) :: ng3,kv3(3,n3d),nstr(n3d)
      INTEGER, INTENT (OUT) :: kxc1_fft,kxc2_fft,kxc3_fft
      INTEGER, INTENT (OUT) :: nxc3_fft,kmxxc_fft
      REAL    bmat(3,3),gmax,gmaxxc,rkmax,sk3(n3d)
c
c---> local variables
c
      INTEGER ksfft,mxc1,mxc2,mxc3,istr,iofile
      REAL    arltv1,arltv2,arltv3,gmxxc_new
c
c---> external subroutines ..
c
      EXTERNAL boxdim
c
c---> intrinsic functions
c
      INTRINSIC int,min
c
c------->          background
c
c        determine the limits  for the xc energy/pot fft-box.
c        since the xc functional has a non-linear dependence 
c        of the charge density in many cases a gmax to determine
c        the xc-potential or energy of gmax=2*rkmax=gmaxq is 
c        in many cases not sufficent. Therefore, we use a larger
c        gmax called gmaxxc to evaluate the xc-potential and energy
c        in the interstitial region. This leads to a larger xc-fft-box.
c        the box dimension is  determined by :
c        m(i) >= gmaxxc*a(i)/(2*pi), where a(i) is the
c        magnitude of the ith basis vector. (this is only true
c        for <a(i)|a(j)> = 0 for i .ne. j., otherwise see boxdim)
c        remember: gmax used through out the flapw code must be
c        larger than or equal  gmax >= gmaxxc >= gmaxp*rkmax.
c
c
c------->          abbreviations
c
c   ksfft=(0,1) : key of selecting fft-prdogram and radix-type
c                      0  program, radix-2 only
c                      1  program, radix-2, radix-3,radix-5
c   iofile      : device number for in and output
c   rkmax       : cut-off for |g+k|
c   gmax        : actually used largest g-vector used througout 
c               : the calculation
c   gmaxxc      : cut-off wavevector for fast fourier transform of 
c                 xc-potential and energy
c                 gmaxxc usually gmaxxc >= 2*rkmax
c   kxc1,2,3d   : dimensions of the xc-pot/energy fft-box
c                 ( in positive fft domain )
c   kxc1,2,3_fft: actual size of the xc-pot/energy fft-box
c                 ( in positive fft domain )
c   nxc3_fft     : number of stars in the xc-pot/energy sphere
c   kmxxc_fft    : number of g-vectors in the xc-pot/energy sphere
c   arltv(i)    : length of reciprical lattice vector along
c                 direction (i)
c
      WRITE (6,'('' gmaxxc should be: 2*kmax <= gmaxxc <= gmax '')')
      IF ( abs( gmaxxc - gmax ) .le. 10.0**(-6) ) THEN
        WRITE (6,'('' concerning memory, you may want to choose'',
     +              '' a smaller value for gmax'')')
      END IF
      IF ( gmaxxc .LE. 10.0**(-6) ) THEN
         WRITE (6,'(" gmaxxc=0 : gmaxxc=gmax choosen as default",
     +              " value")')
         WRITE (6,'(" concerning memory, you may want to choose",
     +              " a smaller value for gmax")')
         gmaxxc=gmax
      END IF
      IF ( gmaxxc .LE. 2*rkmax ) THEN
         WRITE (6,'('' concerning accuracy and total energy'',
     +              '' convergence, you may want'',/,
     +              '' to choose a larger gmaxxc '')')
      END IF
      write (6,'('' gmaxxc ='',f10.6)') gmaxxc
C
C---> Determine dimensions of fft-box of size mxc1, mxc2, mxc3,
c     for which |G(mxc1,mxc2,mxc3)| < Gmaxxc
C
      CALL boxdim(
     >            bmat,
     <            arltv1,arltv2,arltv3)
c
      mxc1 = int( gmaxxc/arltv1 ) + 1
      mxc2 = int( gmaxxc/arltv2 ) + 1
      mxc3 = int( gmaxxc/arltv3 ) + 1
c
c----> fft done in positive domain
c      (remember number of grid points not 2*m+1 because of cyclic
c       boundary condition)
c
      mxc1 = mxc1 + mxc1
      mxc2 = mxc2 + mxc2
      mxc3 = mxc3 + mxc3
c
c---> fft's are usually fastest for low primes
c     (restrict kqid to: kwid=  (2**p) * (3**q) * (5**r)
c
      iofile = 6
      ksfft  = 1
      kxc1_fft = ifft235 (iofile,ksfft,mxc1,2.0)
      kxc2_fft = ifft235 (iofile,ksfft,mxc2,2.0)
      kxc3_fft = ifft235 (iofile,ksfft,mxc3,2.0)
c
c---> for mxc = 2**p, fft is very fast. if mq very close to 2**p
c     choose this, even 2**p < mxc . therefore:
c
      gmxxc_new = min( real( kxc1_fft )*arltv1, 
     +                 real( kxc2_fft )*arltv2, real( kxc3_fft )*arltv3)
      gmxxc_new = 0.5 * gmxxc_new

      IF (gmxxc_new.LT.gmaxxc) THEN
         WRITE (6,'('' gmaxxc recalculated '')')
         WRITE (6,2100) gmaxxc, gmxxc_new, gmxxc_new*gmxxc_new
         WRITE (16,'('' gmaxxc recalculated '')')
         WRITE (16,2100) gmaxxc, gmxxc_new, gmxxc_new*gmxxc_new
         gmaxxc = gmxxc_new
      ENDIF
c
c-----> gmax => gmaxxc
c       (otherwise too few elements in arrays defined in strng1)
c
      IF (gmxxc_new.GT.gmax) THEN
         WRITE (6,'('' gmax must be at least gmxxc_new'')')
         WRITE (6,'('' increase gmax , or reduce gmaxxc'')')
         WRITE (6,'('' gmxxc_new ='',f10.3,''  gmax ='',f10.3)') 
     +                gmxxc_new, gmax
ccc         STOP 'prp_xcfft: gmxxc_new.gt.gmax'
      ENDIF
c
c------> check dimensions of fft chargedensity box used in pwden.f
c
       IF (kxc1_fft.GT.kxc1d .OR. kxc2_fft.gt.kxc2d .OR. 
     +                            kxc3_fft.gt.kxc3d) THEN
          WRITE (6,'('' box dim. for fft too small'')')
          WRITE (6,2110) kxc1_fft,kxc2_fft,kxc3_fft,kxc1d,kxc2d,kxc3d
          WRITE(16,'('' box dim. for fft too small'')')
          WRITE(16,2110) kxc1_fft,kxc2_fft,kxc3_fft,kxc1d,kxc2d,kxc3d
          STOP 'prp_xcfft: mxc[1,2,3]d.gt.kxc[1,2,3]d '
       ENDIF
 2110  FORMAT (' kxc1_fft,kxc2_fft,kxc3_fft,kxc1d,kxc2d,kxc3d ',6i5)
c
c-----> how many stars are in charge density sphere?
c       assume stars are ordered according to length
c
      nxc3_fft = 0
      DO istr = 1 , ng3
         IF ( sk3(istr).LE.gmaxxc ) THEN
            nxc3_fft = istr
         ENDIf
      ENDDO
c
      IF ( nxc3_fft.EQ.0 ) THEN
         WRITE (6,'('' presumably ng3 too small '')')
         WRITE (6,'('' sk3max, gmaxxc '', 2f10.6)')
     +                sk3(ng3),gmaxxc
         WRITE(16,'('' presumably ng3 too small '')')
         WRITE(16,'('' sk3max, gmaxxc '', 2f10.6)')
     +                 sk3(ng3),gmaxxc
         STOP 'prp_xcfft: nxc3_fft.eq.0'
      ENDIF
c
      IF ( nxc3_fft.GT.n3d ) THEN
         WRITE (6,'('' nxc3_fft > n3d '')')
         WRITE (6,'('' nxc3_fft, n3d '',2i10)') nxc3_fft, n3d
         WRITE (16,'('' nxc3_fft > n3d '')')
         WRITE (16,'('' nxc3_fft, n3d '',2i10)') nxc3_fft, n3d
         STOP 'prp_xcfft: nxc3_fft.gt.n3d '
      ENDIF
c
c-----> check that all nxc3_fft stars fit into the xc-pot/energy fft-box
c
      DO istr = 1 , nxc3_fft
        IF ( ( 2*kv3(1,istr).gt.kxc1_fft ) .OR.
     +       ( 2*kv3(2,istr).gt.kxc2_fft ) .OR.
     +       ( 2*kv3(3,istr).gt.kxc3_fft ) ) THEN
          WRITE(6,'('' not all nxc3_fft stars in xc-pot/eng fft box'')')
          WRITE(6,'('' inconsistency in def.s see also strgn1'')')
          WRITE(6,'('' kxc1_fft,kxc2_fft,kxc3_fft,kv1,kv2,kv3 '',6i5)')
     +                 kxc1_fft,kxc2_fft,kxc3_fft,2*kv3(1,istr),
     +                 2*kv3(2,istr),2*kv3(3,istr)
          write(16,'('' not all nxc3_fft stars in xc-pot/eng fft box''
     +              )')
          WRITE(16,'('' inconsistency in def.s see also strgn1'')')
          WRITE(16,'('' kxc1_fft,kxc2_fft,kxc3_fft,kv1,kv2,kv3 '',6i5)')
     +                  kxc1_fft,kxc2_fft,kxc3_fft,2*kv3(1,istr),
     +                  2*kv3(2,istr),2*kv3(3,istr)
          STOP 'prp_xcfft: 2*kv3([1,2,3],istr).gt.mxc[1,2,3]d'
        ENDIF
      ENDDO
c
c-----> how many g-vectors belong to these nxc3_fft stars
c
      kmxxc_fft = 0
      DO istr = 1 , nxc3_fft
         kmxxc_fft = kmxxc_fft + nstr(istr)
      ENDDO

      IF ( kmxxc_fft .gt. kxc1d*kxc2d*kxc3d ) then
         WRITE (6,'('' array dimensions in later subroutines too'',
     +             '' small'')')
      ENDIF

 2100 format(/,1x,'old gmaxxc  =',f10.5, '(a.u.)**(-1) ==>  new gmaxxc'
     >      ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new e_cut(xc)   =',
     >            f10.5,' ry')

      END SUBROUTINE prp_xcfft
      END MODULE m_prpxcfft
