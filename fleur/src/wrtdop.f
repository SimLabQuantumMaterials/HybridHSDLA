      MODULE m_wrtdop
c     ****************************************************
c     write formatted density or potential onto unit 'nu'
c     e. wimmer   march 1985
c     ****************************************************
      CONTAINS
      SUBROUTINE wrtdop(
     >                  jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >                  jspins,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     >                  invs,invs2,film,delz,z1,dx,rmt,zatom,
     >                  nlh,jri,ntypsd,ntypsy,namat,nu,
     >                  iop,dop,it,fr,fpw,fz,fzxy,name)
c
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nu,ntypsd
      INTEGER, INTENT (IN) :: n3d,n2d,nq3,nq2
      INTEGER, INTENT (IN) :: nmzxyd,nmzd,nmzxy,nmz
      INTEGER, INTENT (IN) :: jmtd,nlhd,ntypd,ntype,natd
      INTEGER, INTENT (IN) :: it,nvac,jspd,jspins
      REAL,    INTENT (IN) :: delz,z1
      LOGICAL, INTENT (IN) :: invs,invs2,film
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntypd)
      COMPLEX, INTENT (IN):: fpw(n3d,jspd),fzxy(nmzxyd,n2d-1,2,jspd)
      REAL,    INTENT (IN):: fr(jmtd,0:nlhd,ntypd,jspd),fz(nmzd,2,jspd)
      REAL,    INTENT (IN):: dx(ntypd),rmt(ntypd),zatom(ntypd)
      CHARACTER*8,INTENT (IN):: dop,iop,name(10)
      CHARACTER*2,INTENT (IN):: namat(0:103)
C     ..
C     .. Local Scalars ..
      INTEGER i,ivac,izn,jsp,k,lh,n,na
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
      WRITE (nu) name
      WRITE (6,FMT=8000) name
 8000 FORMAT (' wrtdop title:',10a8)
      WRITE (nu) iop,dop,it
      DO 50 jsp = 1,jspins
         WRITE (nu) jsp
         WRITE (nu) ntype
         na = 1
         DO 20 n = 1,ntype
            izn = zatom(n) + 0.01
            WRITE (nu) namat(izn),n,jri(n),rmt(n),dx(n)
            WRITE (nu) ntypsy(na),nlh(ntypsy(na))
            DO 10 lh = 0,nlh(ntypsy(na))
               WRITE (nu) lh
               WRITE (nu) (fr(i,lh,n,jsp),i=1,jri(n))
   10       CONTINUE
            na = na + neq (n)
   20    CONTINUE
         WRITE (nu) nq3
         IF (invs) THEN
            WRITE (nu) (real(fpw(k,jsp)),k=1,nq3)
         ELSE
            WRITE (nu) (fpw(k,jsp),k=1,nq3)
         END IF
         IF (film) THEN
            DO 40 ivac = 1,nvac
               WRITE (nu) ivac
               WRITE (nu) nmz,z1,delz
               WRITE (nu) (fz(i,ivac,jsp),i=1,nmz)
               WRITE (nu) nq2,nmzxy
               DO 30 k = 2,nq2
                  IF (invs2) THEN
                     WRITE (nu) (real(fzxy(i,k-1,ivac,jsp)),i=1,nmzxy)
                  ELSE
                     WRITE (nu) (fzxy(i,k-1,ivac,jsp),i=1,nmzxy)
                  END IF
   30          CONTINUE
   40       CONTINUE
         END IF
   50 CONTINUE
c
      RETURN
      END SUBROUTINE wrtdop
      END MODULE m_wrtdop
      
