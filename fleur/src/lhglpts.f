      MODULE m_lhglpts
c     **********************************************************
c     calculates lattice harmonics on the gauss-legendre angular 
c     mesh - r.pentcheva Feb'96
c     **********************************************************
      CONTAINS
      SUBROUTINE lhglpts(
     >                   memd,nlhd,ntypsd,lmaxd,nspd,
     >                   rx,nsp,
     >                   clnu,nmem,mlh,nlh,llh,nsymt,
     <                   ylh)
C
      USE m_ylm
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments .. 
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,lmaxd,nspd
      INTEGER, INTENT (IN) :: nsp,nsymt
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: rx(3,nspd)
      REAL,    INTENT (OUT):: ylh(nspd,0:nlhd,ntypsd)
C     ..
C     .. Local Scalars ..
      REAL s
      INTEGER k,lh,mem,nd,ll1,lm
C     ..
C     .. Local Arrays ..
      COMPLEX ylm( (lmaxd+1)**2 )
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
      DO 40 nd = 1,nsymt
         DO 30 k = 1,nsp

            CALL ylm4(
     >                lmaxd,rx(1,k),
     <                ylm)

            DO lh = 0,nlh(nd)
               s = 0
               ll1 = llh(lh,nd) * ( llh(lh,nd) + 1 ) + 1
               DO mem = 1,nmem(lh,nd)
                  lm = ll1 + mlh(mem,lh,nd)
                  s = s + real( clnu(mem,lh,nd) * ylm(lm) )
               ENDDO
               ylh(k,lh,nd) = s
            ENDDO

   30    CONTINUE
   40 CONTINUE
c
      RETURN
      END SUBROUTINE lhglpts
      END MODULE m_lhglpts
