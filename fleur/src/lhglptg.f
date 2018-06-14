      MODULE m_lhglptg
c.....------------------------------------------------------------------
c     calculates lattice harmonics and their gradients on the
c       gauss-legendre angular mesh - r.p. and t.a.
c     for gradient. t.a. 1996.
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE lhglptg(
     >                   nspd,nlhd,ntypsd,lmaxd,memd,
     >                   rx,nsp,igrd,nsymt,
     >                   clnu,nmem,mlh,nlh,llh,
     <                   ylh,thet,ylht1,ylht2,ylhf1,ylhf2,ylhtf)
c
      USE m_polangle
      USE m_ylm
      USE m_dylm

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nspd,nlhd,ntypsd,lmaxd,memd
      INTEGER, INTENT (IN) :: nsp,igrd,nsymt
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: rx(3,nspd)
      REAL,    INTENT (OUT):: ylh(nspd,0:nlhd,ntypsd),thet(nspd)
      REAL,    INTENT (OUT):: ylht1(nspd,0:nlhd,ntypsd)
      REAL,    INTENT (OUT):: ylht2(nspd,0:nlhd,ntypsd)
      REAL,    INTENT (OUT):: ylhtf(nspd,0:nlhd,ntypsd)
      REAL,    INTENT (OUT):: ylhf1(nspd,0:nlhd,ntypsd)
      REAL,    INTENT (OUT):: ylhf2(nspd,0:nlhd,ntypsd)
C     ..
C     .. Local Scalars ..
      REAL s,st1,st2,sf1,sf2,stf,phi
      INTEGER k,lh,mem,nd,lm,ll1
C     ..
C     .. Local Arrays ..
      COMPLEX ylm( (lmaxd+1)**2 )
      COMPLEX dylmt1( (lmaxd+1)**2 ), dylmt2( (lmaxd+1)**2 )
      COMPLEX dylmf1( (lmaxd+1)**2 ), dylmf2( (lmaxd+1)**2 )
      COMPLEX dylmtf( (lmaxd+1)**2 )
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
c.....------------------------------------------------------------------
c     ..
      DO 40 nd = 1,nsymt

          DO 30 k = 1,nsp

              CALL ylm4(
     >                   lmaxd,rx(1,k),
     <                   ylm)
              CALL pol_angle(
     >                       rx(1,k),rx(2,k),rx(3,k),
     <                       thet(k),phi)

              IF (igrd.GT.0) THEN
                CALL dylm3(
     >                     lmaxd,lmaxd,rx(1,k),ylm,
     <                     dylmt1,dylmt2,dylmf1,dylmf2,dylmtf)
              ENDIF

              DO 20 lh = 0,nlh(nd)
                  s   = 0
                  st1 = 0
                  st2 = 0
                  sf1 = 0
                  sf2 = 0
                  stf = 0
                  ll1 = llh(lh,nd) * ( llh(lh,nd) + 1 ) + 1

                  DO mem = 1,nmem(lh,nd)
                    lm = ll1 + mlh(mem,lh,nd)
                    s = s + real( clnu(mem,lh,nd) * ylm(lm) )
                  ENDDO

                  ylh(k,lh,nd) = s

                  IF (igrd.GT.0) THEN

                    DO mem = 1,nmem(lh,nd)
                      lm = ll1 + mlh(mem,lh,nd)
                      s   = s   + real( clnu(mem,lh,nd)* ylm(lm) )
                      st1 = st1 + real( clnu(mem,lh,nd)*dylmt1(lm) )
                      st2 = st2 + real( clnu(mem,lh,nd)*dylmt2(lm) )
                      sf1 = sf1 + real( clnu(mem,lh,nd)*dylmf1(lm) )
                      sf2 = sf2 + real( clnu(mem,lh,nd)*dylmf2(lm) )
                      stf = stf + real( clnu(mem,lh,nd)*dylmtf(lm) )
                    ENDDO

                    ylht1(k,lh,nd) = st1
                    ylht2(k,lh,nd) = st2
                    ylhf1(k,lh,nd) = sf1
                    ylhf2(k,lh,nd) = sf2
                    ylhtf(k,lh,nd) = stf

                  ENDIF

   20         ENDDO
   30     ENDDO
   40 ENDDO

      RETURN
      END SUBROUTINE lhglptg
      END MODULE m_lhglptg
