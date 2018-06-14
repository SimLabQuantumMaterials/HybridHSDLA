      MODULE m_od_strgn1
      CONTAINS
      SUBROUTINE od_strgn1(
     >     igrd,bmat,zrfs,invs,invs2,odd,
     <     ig,kv,nstr,
     <     igfft,pgfft,pgftx,pgftxx,pgftxy,pgfty,pgftyy)

c**********************************************************
c     one-dimensional stars generator
c                                 Y.Mokrousov, 2002-2003
c***********************************************************
      USE m_constants, ONLY : pimach
      USE m_od_types, ONLY : od_dim

      IMPLICIT NONE

      TYPE (od_dim), INTENT (IN) :: odd

      INTEGER, INTENT (IN) :: igrd
      REAL,    INTENT (IN) :: bmat(3,3)
      LOGICAL, INTENT (IN) :: zrfs,invs,invs2

      INTEGER, INTENT (OUT) :: ig(-odd%k3:odd%k3,-odd%M:odd%M)
      INTEGER, INTENT (OUT) :: kv(2,odd%n2d),nstr(odd%n2d)
      INTEGER, INTENT (OUT) :: igfft(0:odd%nn2d-1,2)
      REAL,    INTENT (OUT) :: pgfft(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgftx(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgftxx(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgftxy(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgfty(0:odd%nn2d-1)
      REAL,    INTENT (OUT) :: pgftyy(0:odd%nn2d-1)

      INTEGER nfftx_1,nffty_1,kfx_1,kfy_1,kfft_1
      REAL gfx_1,gfy_1
      INTEGER i,m,m1,z,z1

c     odl%nn2d = odi%nn2d
c     odg%nn2d = odi%nn2d

c     ALLOCATE ( odi%kv(1:2,odi%n2d),
c    &           odi%ig(-odi%k3:odi%k3,-odi%M:odi%M),
c    &           odi%nst2(1:odi%n2d),
c    &           odl%igf(0:odl%nn2d-1,1:2),
c    &           odl%pgf(0:odl%nn2d-1) )
      igfft(0:odd%nn2d-1,1:2) = 0
      pgfft(0:odd%nn2d-1) = 0.

      IF (igrd.NE.0) THEN
c         ALLOCATE ( odg%pgfx(0:odg%nn2d-1),
c     &              odg%pgfy(0:odg%nn2d-1),
c     &              odg%pgfxx(0:odg%nn2d-1),
c     &              odg%pgfxy(0:odg%nn2d-1),
c     &              odg%pgfyy(0:odg%nn2d-1) )
         pgftx(0:odd%nn2d-1) = 0.
         pgfty(0:odd%nn2d-1) = 0.
         pgftxx(0:odd%nn2d-1) = 0.
         pgftxy(0:odd%nn2d-1) = 0.
         pgftyy(0:odd%nn2d-1) = 0.
      END IF

c      odi%nq2 = odi%n2d
c      odd%kimax2 = odd%n2d - 1
      nstr(1:odd%n2d) = 0

c---> generating mapping arrays       

      i = 0

      DO 11 z1 = 0,2*odd%k3
         DO 12 m1 = 0,2*odd%M
            IF (m1.LE.odd%M) m = m1
            IF (m1.GT.odd%M) m = odd%M - m1
            IF (z1.LE.odd%k3) z = z1
            IF (z1.GT.odd%k3) z = odd%k3-z1
            IF (odd%chi.EQ.1) THEN
               IF (zrfs) THEN
                  IF (MOD(m,odd%rot).EQ.0) THEN
                     IF (z.GE.0) THEN
                        i = i+1
                        ig(z,m) = i
                        kv(1,i) = z
                        kv(2,i) = m
                        IF (z.EQ.0) THEN
                           nstr(i) = 1
                        ELSE
                           nstr(i) = 2
                        END IF
                     ELSE
                        ig(z,m) = ig(-z,m)
                     END IF
                  ELSE
                     ig(z,m) = 0
                  END IF
               ELSE IF (.NOT.zrfs) THEN
                  IF (MOD(m,odd%rot).EQ.0) THEN
                     i = i+1
                     ig(z,m) = i
                     kv(1,i) = z
                     kv(2,i) = m
                     nstr(i) = 1
                  ELSE
                     ig(z,m) = 0
                  END IF
               END IF
            ELSE
               IF (MOD(m+(odd%rot)*z,odd%chi).EQ.0) THEN
                  i = i+1
                  ig(z,m) = i
                  kv(1,i) = z
                  kv(2,i) = m
                  nstr(i) = 1
               ELSE
                  ig(z,m) = 0 
               END IF
            END IF
 12      CONTINUE
 11   CONTINUE

c---> preparations for 2dim vacuum fft
c---> at the moment we have no symmetries

      nfftx_1 = 3*odd%k3
      nffty_1 = 3*odd%M
      DO 2001 i = 1,odd%nq2
         kfx_1 = kv(1,i)
         kfy_1 = kv(2,i)
         IF (kfx_1.LT.0) kfx_1 = kfx_1 + nfftx_1
         IF (kfy_1.LT.0) kfy_1 = kfy_1 + nffty_1
         kfft_1 = kfx_1 + kfy_1*nfftx_1
         igfft(i-1,1) = i
         igfft(i-1,2) = kfft_1
         pgfft(i-1) = 1.
 2001 END DO

      IF (igrd.NE.0) THEN
         DO 2002 i = 1,odd%nq2
            kfx_1 = kv(1,i)
            kfy_1 = kv(2,i)
            gfx_1 = bmat(3,3)*kfx_1
            gfy_1 = kfy_1
            pgftx(i-1)  = gfx_1
            pgfty(i-1)  = gfy_1
            pgftxx(i-1) = gfx_1*gfx_1
            pgftyy(i-1) = gfy_1*gfy_1
            pgftxy(i-1) = gfx_1*gfy_1
 2002    CONTINUE
      END IF

c     odi%kimax2 = odi%nq2 - 1
      
      RETURN
      END SUBROUTINE od_strgn1
      END MODULE m_od_strgn1
 

