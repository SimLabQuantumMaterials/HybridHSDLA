      MODULE m_tlo
c***********************************************************************
c     sets up the extra t-matrix elements due to the local orbitals.
c     only non=zero elements are calculated
c
c     p.kurz jul. 1996
c***********************************************************************
      CONTAINS
      SUBROUTINE tlo(
     >               lmaxd,ntypd,jspd,nwdd,jmtd,nlod,llod,nlhd,ntypsd,
     >               lmd,loplod,memd,nlh,llh,nmem,mlh,clnu,ntypsy,
     >               jsp,ntyp,lmx,jri,r0,dx,el,ello,
     >               lh0,secvar,vr,natd,nat,l_dulo,ulo_der,
     >               flo,f,g,nlo,llo,lo1l,uulon,dulon,uloulopn,
     >               uuilon,duilon,ulouilopn,
     <               tuulo,tdulo,tuloulo)
c
c*************** ABBREVIATIONS *****************************************
c tuulo      : t-matrix element of the lo and the apw radial fuction
c tdulo      : t-matrix element of the lo and the energy derivativ of 
c              the apw radial fuction
c tuloulo    : t-matrix element of two los
cc***********************************************************************
c
      USE m_intgr, ONLY : intgr3  
      USE m_gaunt, ONLY: gaunt1
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,jspd,nwdd,jmtd,nlod,llod
      INTEGER, INTENT (IN) :: lmd,loplod,nlhd,ntypsd,memd,natd
      INTEGER, INTENT (IN) :: jri,jsp,lmx,ntyp,nat,lh0
      REAL,    INTENT (IN) :: dx
      LOGICAL, INTENT (IN) :: secvar
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),ntypsy(natd) 
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),lo1l(0:llod,ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: el(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,jspd),vr(jmtd,0:nlhd)
      REAL,    INTENT (IN) :: f(jmtd,2,0:lmaxd),g(jmtd,2,0:lmaxd)
      REAL,    INTENT (IN) :: flo(jmtd,2,nlod),uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN) :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      REAL,    INTENT (IN) :: uuilon(nlod,ntypd),duilon(nlod,ntypd)
      REAL,    INTENT (IN) :: ulouilopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN) :: r0(jmtd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (OUT) :: tdulo(0:lmd,-llod:llod,*)
      COMPLEX, INTENT (OUT) :: tuulo(0:lmd,-llod:llod,*)
      COMPLEX, INTENT (OUT) :: tuloulo(-llod:llod,-llod:llod,*)
C     ..
C     .. Local Scalars ..
      COMPLEX ci,cil
      INTEGER i,l,lh,lm,lmax,lmin,lmp,lo,lop,loplo,lp,lpmax,lpmax0,
     +        lpmin,lpmin0,lpp,m,mem,mp,mpp
C     ..
C     .. Local Arrays ..
      REAL x(jmtd),ulovulo(loplod,lh0:nlhd)
      REAL uvulo(nlod,0:lmaxd,lh0:nlhd),dvulo(nlod,0:lmaxd,lh0:nlhd)
C     ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,max,min,mod
C     ..
      ci = cmplx(0.0,1.0)
      DO lo = 1,nlo(ntyp)
         l = llo(lo,ntyp)
         DO lp = 0,lmx
            DO lh = lh0,nlh(ntypsy(nat))
               lpp = llh(lh,ntypsy(nat))
               lmin = abs(lp-l)
               lmin = lp - l
               lmax = lp + l
               IF ((mod(l+lp+lpp,2).EQ.1) .OR. (lpp.LT.lmin) .OR.
     +             (lpp.GT.lmax)) THEN
                  uvulo(lo,lp,lh) = 0.0
                  dvulo(lo,lp,lh) = 0.0
               ELSE
                  DO i = 1,jri
                     x(i) = (f(i,1,lp)*flo(i,1,lo)+
     +                      f(i,2,lp)*flo(i,2,lo))*vr(i,lh)
                  END DO
                  CALL intgr3(x,r0,dx,jri,uvulo(lo,lp,lh))
                  DO i = 1,jri
                     x(i) = (g(i,1,lp)*flo(i,1,lo)+
     +                      g(i,2,lp)*flo(i,2,lo))*vr(i,lh)
                  END DO
                  CALL intgr3(x,r0,dx,jri,dvulo(lo,lp,lh))
               END IF
            END DO
         END DO
      END DO
      loplo = 0
      DO lop = 1,nlo(ntyp)
         lp = llo(lop,ntyp)
         DO lo = 1,lop
            l = llo(lo,ntyp)
            loplo = loplo + 1
            IF (loplo.GT.loplod) STOP 'tlo: loplo > loplod!!!'
            DO lh = lh0,nlh(ntypsy(nat))
               lpp = llh(lh,ntypsy(nat))
               lmin = abs(lp - l)
               lmax = lp + l
               IF ((mod(l+lp+lpp,2).EQ.1) .OR. (lpp.LT.lmin) .OR.
     +             (lpp.GT.lmax)) THEN
                  ulovulo(loplo,lh) = 0.0
               ELSE
                  DO i = 1,jri
                     x(i) = (flo(i,1,lop)*flo(i,1,lo)+
     +                      flo(i,2,lop)*flo(i,2,lo))*vr(i,lh)
                  END DO
                  CALL intgr3(x,r0,dx,jri,ulovulo(loplo,lh))
               END IF
            END DO
         END DO
      END DO
c---> generate the different t matrices
c---> but first initialize them ( done in eigen )
!     
c---> generate the t-matrices. for optimal performance consider only
c---> those combinations of l,l',l'',m,m',m'' that satisfy the three
c---> conditions for non-zero gaunt-coeff. i.e.
c---> |l - l''| <= l' <= l + l'' (triangular condition)
c---> m' = m + m'' and l + l' + l'' even
c---> loop over the local orbitals
      DO lo = 1,nlo(ntyp)
         l = llo(lo,ntyp)
         DO m = -l,l
c--->       loop over the lattice harmonics
            DO lh = lh0,nlh(ntypsy(nat))
               lpp = llh(lh,ntypsy(nat))
               lpmin0 = abs(l-lpp)
               lpmax0 = l + lpp
c--->          check that lpmax is smaller than the max l of the
c--->          wavefunction expansion at this atom
               lpmax = min(lpmax0,lmx)
c--->          make sure that l + l'' + lpmax is even
               lpmax = lpmax - mod(l+lpp+lpmax,2)
               DO mem = 1,nmem(lh,ntypsy(nat))
                  mpp = mlh(mem,lh,ntypsy(nat))
                  mp = m + mpp
                  lpmin = max(lpmin0,abs(mp))
c--->             make sure that l + l'' + lpmin is even
                  lpmin = lpmin + mod(abs(lpmax-lpmin),2)
c--->             loop over l'
                  DO lp = lpmin,lpmax,2
                     lmp = lp* (lp+1) + mp
                     cil = ((ci** (l-lp))*clnu(mem,lh,ntypsy(nat)))*
     +                     gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)
                     tuulo(lmp,m,lo) = tuulo(lmp,m,lo) +
     +                                 cil*uvulo(lo,lp,lh)
                     tdulo(lmp,m,lo) = tdulo(lmp,m,lo) +
     +                                 cil*dvulo(lo,lp,lh)
                  END DO
               END DO
            END DO
         END DO
      END DO
c---> generate the t-matrix including two local orbitals for lo' >= lo
c---> loop over lo'
      loplo = 0
      DO lop = 1,nlo(ntyp)
         lp = llo(lop,ntyp)
         DO mp = -lp,lp
c--->       loop over the lattice harmonics
            DO lh = lh0,nlh(ntypsy(nat))
               lpp = llh(lh,ntypsy(nat))
               DO mem = 1,nmem(lh,ntypsy(nat))
                  mpp = mlh(mem,lh,ntypsy(nat))
                  m = mp - mpp
c--->             loop over lo
                  DO lo = 1,lop
                     l = llo(lo,ntyp)
                     loplo = (lop-1)*lop/2 + lo
                     IF ((abs(l-lpp).LE.lp) .AND. (lp.LE. (l+lpp)) .AND.
     +                   (mod(l+lp+lpp,2).EQ.0) .AND.
     +                   (abs(m).LE.l)) THEN
                        cil = ((ci** (l-lp))*clnu(mem,lh,ntypsy(nat)))*
     +                        gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)
                        tuloulo(mp,m,loplo) = tuloulo(mp,m,loplo) +
     +                                        cil*ulovulo(loplo,lh)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
c---> add the diagonal terms from the muffin-tin hamiltonian. these
c---> terms have to be made hermitian. if second variation is switched
c---> on, the t-matrices contain only the contributions from the
c---> non-spherical hamiltonian.
      IF (.NOT.secvar) THEN
         DO lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            DO m = -l,l
               lm = l* (l+1) + m
               tuulo(lm,m,lo) = tuulo(lm,m,lo) + 0.5 * uulon(lo,ntyp) *
     +                           ( el(l,ntyp,jsp,1)+ello(lo,ntyp,jsp) )
               tdulo(lm,m,lo) = tdulo(lm,m,lo) + 0.5 * dulon(lo,ntyp) *
     +                           ( el(l,ntyp,jsp,1)+ello(lo,ntyp,jsp) )
     +                           + 0.5 * uulon(lo,ntyp)
               IF (ulo_der(lo,ntyp).ge.1) THEN
                 tuulo(lm,m,lo) = tuulo(lm,m,lo) + 0.5 * uuilon(lo,ntyp)
                 tdulo(lm,m,lo) = tdulo(lm,m,lo) + 0.5 * duilon(lo,ntyp)
               ENDIF
c+apw_lo
               IF (l_dulo(lo,ntyp)) THEN         
                 tuulo(lm,m,lo) = tuulo(lm,m,lo) + 0.5
                 tdulo(lm,m,lo) = 0.0
               ENDIF
c+apw_lo
c--->          the 1 in el(l,ntyp,jsp,1) is the window index. the
c--->          programm is not intended to do multiple window
c--->          calculations and local orbitals at the same time.
            END DO
         END DO
         DO lop = 1,nlo(ntyp)
            lp = llo(lop,ntyp)
            DO lo = lo1l(lp,ntyp),lop
               loplo = (lop-1)*lop/2 + lo
               DO m = -lp,lp
                  tuloulo(m,m,loplo) = tuloulo(m,m,loplo) +
     +                                 0.5* (ello(lop,ntyp,jsp)+
     +                                 ello(lo,ntyp,jsp))*
     +                                 uloulopn(lop,lo,ntyp) +
     +                                 0.5* (
     +                                 ulouilopn(lop,lo,ntyp) +
     +                                 ulouilopn(lo,lop,ntyp))
               END DO
            END DO
         END DO
      END IF

      END SUBROUTINE tlo
      END MODULE m_tlo
