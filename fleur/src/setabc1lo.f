      MODULE m_setabc1lo
c*********************************************************************
c calculates the (lower case) a, b and c coefficients for the local
c orbitals. The radial function of the local orbital is a linear 
c combination of the apw radial function and its derivative and the
c extra radial funtion (a*u + b*udot + c*ulo). This function is zero
c and has zero derivative at the muffin tin boundary.
c Philipp Kurz 99/04
c*********************************************************************
      CONTAINS
      SUBROUTINE setabc1lo(
     >                     ntypd,lmaxd,nlod,
     >                     ntyp,us,dus,uds,duds,ddn,
     >                     nlo,llo,l_dulo,ulos,dulos,uulon,dulon,
     <                     alo1,blo1,clo1)
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: ntypd,lmaxd,nlod
      INTEGER, INTENT (IN)  :: ntyp
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN)  :: us(0:lmaxd,ntypd),dus(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: uds(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: ulos(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (IN)  :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      REAL,    INTENT (IN)  :: ddn(0:lmaxd,ntypd)
      LOGICAL, INTENT (IN)  :: l_dulo(nlod,ntypd)
      REAL,    INTENT (OUT) :: alo1(nlod),blo1(nlod),clo1(nlod)
C     ..
C     .. Local Scalars ..
      REAL ka,kb,ws
      INTEGER l,lo
      LOGICAL apw_at
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
c look, whether 'ntyp' is a APW atom; then set apw_at=.true.
c
      apw_at = .false.
      DO lo = 1,nlo(ntyp)
         IF (l_dulo(lo,ntyp)) apw_at = .true.
      ENDDO
      DO lo = 1,nlo(ntyp)
         l = llo(lo,ntyp)
         IF (apw_at) THEN
           IF (l_dulo(lo,ntyp)) THEN
c udot lo
             ka = sqrt( 1+(us(l,ntyp)/uds(l,ntyp))**2 * ddn(l,ntyp))
             alo1(lo)=1.00 / ka
             blo1(lo)=-us(l,ntyp)/ (uds(l,ntyp) * ka )
             clo1(lo)=0.00
           ELSE
c u2 lo
             alo1(lo)=1.00
             blo1(lo)=0.00
             clo1(lo)=-us(l,ntyp)/ulos(lo,ntyp)           
           ENDIF
         ELSE
           ws = uds(l,ntyp)*dus(l,ntyp) - us(l,ntyp)*duds(l,ntyp)
           ka = 1.0/ws* (duds(l,ntyp)*ulos(lo,ntyp)-
     +          uds(l,ntyp)*dulos(lo,ntyp))
           kb = 1.0/ws* (us(l,ntyp)*dulos(lo,ntyp)-
     +          dus(l,ntyp)*ulos(lo,ntyp))
           clo1(lo) = 1.0/sqrt(ka**2+ (kb**2)*ddn(l,ntyp)+1.0+
     +                2.0*ka*uulon(lo,ntyp)+2.0*kb*dulon(lo,ntyp))
           alo1(lo) = ka*clo1(lo)
           blo1(lo) = kb*clo1(lo)
         ENDIF 
      END DO

      END SUBROUTINE setabc1lo
      END MODULE m_setabc1lo
c
c         flo = alo1(lo)*us(l,ntyp) + blo1(lo)*uds(l,ntyp) +
c     +         clo1(lo)*ulos(lo,ntyp)
c         dflo = alo1(lo)*dus(l,ntyp) + blo1(lo)*duds(l,ntyp) +
c     +          clo1(lo)*dulos(lo,ntyp)
c         nflo = alo1(lo)**2 + (blo1(lo)**2)*ddn(l,ntyp) + clo1(lo)**2 +
c     +          2*alo1(lo)*clo1(lo)*uulon(lo,ntyp) +
c     +          2*blo1(lo)*clo1(lo)*dulon(lo,ntyp)
