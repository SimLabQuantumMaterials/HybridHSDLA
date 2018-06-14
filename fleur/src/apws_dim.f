      SUBROUTINE apws_dim(
     >                    bkpt,bmat,bbmat,rkmax,l_ss,qss,odd,
     <                    nv,nv2,kq1d,kq2d,kq3d)
c
c*********************************************************************
c     determines the lapw list such that |k+G|<rkmax.
c     bk(i) is the nk k-point given in internal (i.e. b1,b2,b3) units.
c     the lapw list is stored on unit29 in direct access mode.
c     (it is the iteration number and nk0 is the base value
c     for the k-point list, including windows.)
c        m. weinert  1986
c
c     dimensions kq(i)d for charge density FFT added.
c        s. bluegel, JRCAT, Feb. 97
c     subroutine boxdim
c        s.bluegel, IFF, 18.Nov.97^F
c*********************************************************************
      USE m_ifft,     ONLY : ifft235
      USE m_dotir,    ONLY : dotirp
      USE m_od_types, ONLY : od_dim

      IMPLICIT NONE
      INTEGER nv,nv2,kq1d,kq2d,kq3d
      INTEGER j1,j2,j3,mk1,mk2,mk3,iofile,ksfft
      INTEGER ispin,nvh(2),nv2h(2)
      LOGICAL l_ss
      REAL rkmax,bkpt(3),bmat(3,3),bbmat(3,3),qss(3)
      REAL arltv1,arltv2,arltv3,rkm,rk2,r2,s(3),gmaxp
c-odim
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C ..
C ... External Subroutines
      EXTERNAL boxdim
C ..
C ... Intrinsic Functions
      INTRINSIC int
C
C------->          ABBREVIATIONS
C
C   ksfft=(0,1) : KEY OF SELECTING FFT-PRDOGRAM AND RADIX-TYPE
C                      0  PROGRAM, RADIX-2 ONLY
C                      1  PROGRAM, RADIX-2, RADIX-3,RADIX-5
c   iofile      : device number for in and output
c   gmax        : cut-off wavevector for charge density
c   rkmax       : cut-off for |g+k|
c   gmaxp       : gmaxp = gmax/rkmax, ideal: gmaxp=2
c   arltv(i)    : length of reciprical lattice vector along
c                 direction (i)
c
C---> Determine rkmax box of size mk1, mk2, mk3,
c     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
C
      CALL boxdim(
     >            bmat,
     <            arltv1,arltv2,arltv3)

c     (add 1+1 due to integer rounding, strange k_vector in BZ)
      mk1 = int(rkmax/arltv1) + 2
      mk2 = int(rkmax/arltv2) + 2
      mk3 = int(rkmax/arltv3) + 2

      rkm = rkmax
      rk2 = rkm*rkm
c---> obtain vectors
c---> in a spin-spiral calculation different basis sets are used for
c---> the two spin directions, because the cutoff radius is defined
c---> by |G + k +/- qss/2| < rkmax.
      DO ispin = 1,2
         nv = 0
         nv2 = 0
         DO j1 = -mk1,mk1
            s(1) = bkpt(1) + j1 + (2*ispin - 3)/2.0*qss(1)
            DO j2 = -mk2,mk2
               s(2) = bkpt(2) + j2 + (2*ispin - 3)/2.0*qss(2)
c--->          nv2 for films
               s(3) = 0.0
               r2 = dotirp(s,s,bbmat)
               IF (r2.LE.rk2) nv2 = nv2 + 1
               DO j3 = -mk3,mk3
                  s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*qss(3)
                  r2 = dotirp(s,s,bbmat)
                  IF (r2.LE.rk2) THEN
                     nv = nv + 1
                  END IF
               END DO
            END DO
         END DO
c-odim
         IF (odd%d1) THEN
           nv2 = 0
           s(1) = 0.0
           s(2) = 0.0
           DO j3 = -mk3,mk3
              s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*qss(3)
              r2 = dotirp(s,s,bbmat)
              IF (r2.LE.rk2) THEN
                 nv2 = nv2 + 1
              END IF
           END DO
         END IF
c+odim
         nvh(ispin)  = nv
         nv2h(ispin) = nv2
      END DO
      nv  = max(nvh(1),nvh(2))
      nv2 = max(nv2h(1),nv2h(2))

c---> Determine the dimensions kq1d, kq2d, kq3d
c     of the dimension of the charge density fft-box
c     needed for the fast calculation of pw density
c     (add 1 due to integer rounding,
c      factor 2 due to positive domain)
c
      gmaxp = 2.0
      CALL boxdim(bmat,arltv1,arltv2,arltv3)
c
      mk1 = int(gmaxp*rkmax/arltv1) + 1
      mk2 = int(gmaxp*rkmax/arltv2) + 1
      mk3 = int(gmaxp*rkmax/arltv3) + 1

c---> add + 1 in spin spiral calculation, to make sure that all G's are 
c---> still within the FFT-box after being shifted by the spin spiral
c---> q-vector.
      IF (l_ss) THEN
         mk1 = mk1 + 1
         mk2 = mk2 + 1
         mk3 = mk3 + 1
      ENDIF         
c
      kq1d = 2*mk1
      kq2d = 2*mk2
      kq3d = 2*mk3
c
c---> fft's are usually fastest for low primes
c     (restrict kqid to: kqid=  (2**P) * (3**Q) * (5**R)
c
      iofile = 6
      ksfft = 1
c      WRITE (6,*) 'minimum: kq1d,kq2d,kq3d',kq1d,kq2d,kq3d
      kq1d = ifft235(iofile,ksfft,kq1d,gmaxp)
      kq2d = ifft235(iofile,ksfft,kq2d,gmaxp)
      kq3d = ifft235(iofile,ksfft,kq3d,gmaxp)
c      WRITE (6,*) 'ifft235: kq1d,kq2d,kq3d',kq1d,kq2d,kq3d

      RETURN
      END
