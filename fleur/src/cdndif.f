      SUBROUTINE cdndif(
     >                  m1,m2,m3,n1,n2,n3,s1,s2,film,natd,neq,
     >                  n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd,
     >                  nq3,nq2,jri,nlh,ntype,nmz,nmzxy,ntypsy,nvac,
     X                  rhsp,rhpw,rhv0,rhv1)
c
c     ******************************************************
c     forms the difference or the addition of two density
c     sets                       c.l.fu
c     ******************************************************

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,n2d,jmtd,nlhd,ntypd,nmzd,nmzxyd,ntypsd
      INTEGER, INTENT (IN) :: nq3,nq2,ntype,nmz,nmzxy,nvac,natd
      INTEGER, INTENT (IN) :: m1,m2,m3,n1,n2,n3
      REAL,    INTENT (IN) :: s1,s2
      LOGICAL, INTENT (IN) :: film
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: jri(ntypd)
      INTEGER, INTENT (IN)    :: ntypsy(natd),neq(ntypd),nlh(ntypsd)
      COMPLEX, INTENT (INOUT) :: rhpw(n3d,2,2)
      COMPLEX, INTENT (INOUT) :: rhv1(nmzxyd,n2d-1,2,2,2)
      REAL,    INTENT (INOUT) :: rhsp(jmtd,0:nlhd,ntypd,2,2)
      REAL,    INTENT (INOUT) :: rhv0(nmzd,2,2,2)
C     ..
C     .. Local Scalars ..
      INTEGER ivac,j,k,lh,n,na
C     ..
      na = 1
      DO 30 n = 1,ntype
         DO 20 lh = 0,nlh(ntypsy(na))
            DO 10 j = 1,jri(n)
               rhsp(j,lh,n,m1,n1) = s1*rhsp(j,lh,n,m2,n2) +
     +                              s2*rhsp(j,lh,n,m3,n3)
   10       CONTINUE
   20    CONTINUE
         na = na + neq(n)
   30 CONTINUE
      DO 40 k = 1,nq3
         rhpw(k,m1,n1) = s1*rhpw(k,m2,n2) + s2*rhpw(k,m3,n3)
   40 CONTINUE
      IF (film) THEN
         DO 80 ivac = 1,nvac
            DO 50 j = 1,nmz
               rhv0(j,ivac,m1,n1) = s1*rhv0(j,ivac,m2,n2) +
     +                              s2*rhv0(j,ivac,m3,n3)
   50       CONTINUE
            DO 70 j = 1,nmzxy
               DO 60 k = 1,nq2 - 1
                  rhv1(j,k,ivac,m1,n1) = s1*rhv1(j,k,ivac,m2,n2) +
     +                                   s2*rhv1(j,k,ivac,m3,n3)
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
      END IF
c      write(16,*)'in cdndif, num, id=',m1,n1,'s1,s2',s1,s2
c      write(16,'(8d10.4)')(rhsp(i,0,1,m1,n1),i=1,20)
      RETURN
      END
