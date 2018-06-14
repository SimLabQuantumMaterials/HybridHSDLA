      MODULE m_convn
      CONTAINS
      SUBROUTINE convn(
     >                 ncvd,ntype,gmax,rmt,
     <                 ncv)
c
c     ***********************************************************
c     determines the optimum values for the convergence parameter
c     for each atom type using the criterion discussed in
c     m. weinert, j. math. phys. 22, 2433 (1981).  each sphere
c     and l component may have different values.  (psqpw changed
c     to allow this option).
c          m. weinert july 1982
c     ***********************************************************
      IMPLICIT NONE
C     ..
C     .. Scalars Arguments ..
      INTEGER, INTENT (IN)  :: ntype,ncvd
      REAL,    INTENT (IN)  :: gmax
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (OUT) :: ncv(ntype)
      REAL,    INTENT (IN)  :: rmt(ntype)
C     ..
C     .. Local Scalars ..
      REAL sck,z0
      INTEGER i,l,n,n1,nc
C     ..
C     .. Local Arrays ..
      REAL z(17)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC min0
C     ..
C     .. Data statements ..
      DATA z/6.9e0,8.1e0,9.3e0,10.5e0,11.6e0,12.7e0,13.9e0,15.0e0,
     +     16.1e0,17.2e0,18.3e0,19.4e0,20.5e0,21.6e0,22.7e0,23.7e0,
     +     24.8e0/,z0/5.7e0/
C     ..
c--->    read in values of ncv (if ncv(1).le.0, calculate best values)
c      read(5,1000) (ncv(n),n=1,ntype)
c      if(ncv(1).le.0) go to 2
c      n1=ncv(1)
c      do 1 n=2,ntype
c    1 if(ncv(n).le.0) ncv(n)=n1
c      go to 5
c--->    calculate values
c    2 continue
c
      DO 20 n = 1,ntype
         sck = gmax*rmt(n)
         IF (sck.LT.z0) GO TO 60
         DO 10 i = 1,17
            IF (sck.GT.z(i)) GO TO 10
            ncv(n) = i
            GO TO 20
   10    CONTINUE
         n1 = 0.9e0* (sck-z(17))
         ncv(n) = 18 + n1
   20 CONTINUE
c--->    output and make sure ncv(n).le.ncvd
   30 CONTINUE
      WRITE (6,FMT=8010)
      WRITE (16,FMT=8010)
      DO 40 n = 1,ntype
         nc = ncv(n)
         l = nc - 1
         WRITE (6,FMT=8020) n,nc,l
         WRITE (16,FMT=8020) n,nc,l
   40 CONTINUE
      l = ncvd - 1
      WRITE (6,FMT=8030) ncvd,l
      WRITE (16,FMT=8030) ncvd,l
      DO 50 n = 1,ntype
         ncv(n) = min0(ncv(n),ncvd)
   50 CONTINUE
      RETURN
   60 WRITE (6,FMT=8040) n,sck
      WRITE (16,FMT=8040) n,sck
      STOP 'ncv'
 8000 FORMAT (10i5)
 8010 FORMAT (/,/,10x,'convergence parameters for the pseudocharge',
     +       ' density expansion',/,10x,'atom',5x,'parameter',5x,
     +       'max. l to include',/)
 8020 FORMAT (10x,i3,9x,i3,13x,i3)
 8030 FORMAT (10x,'max values allowed: ncvd=',i3,', l=',i3,/)
 8040 FORMAT (/,/,10x,'atom type',i3,' has rkmax=',f6.4,/,10x,
     +       '$$$ stop ncv error')
      END SUBROUTINE convn
      END MODULE m_convn
