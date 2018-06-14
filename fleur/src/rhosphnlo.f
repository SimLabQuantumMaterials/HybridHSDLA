      MODULE m_rhosphnlo
c***********************************************************************
c Add the local orbital contributions to the charge density. The 
c corresponding summation of the pure apw contribuions is done in
c cdnval.
c Philipp Kurz 99/04
c***********************************************************************
      CONTAINS
      SUBROUTINE rhosphnlo(
     >                     nlod,lmaxd,nlhd,jmtd,ntypsd,natd,
     >                     ntypsy,nlh,uloulopn,dulon,uulon,llo,nlo,nat,
     >                     neq,lmax,ello,vr,jri,r0,dx,sfp,
     >                     aclo,bclo,cclo,acnmt,bcnmt,ccnmt,f,g,l_dulo,
     >                     ulo_der,
     X                     rho,qmtllo)

      USE m_constants, ONLY : c_light
      USE m_radsra
      USE m_radsrdn

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nlod,lmaxd,nlhd,jmtd,ntypsd
      INTEGER, INTENT (IN) :: jri,nat,natd,lmax,nlo,neq
      REAL,    INTENT (IN) :: r0,dx,sfp
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod),nlh(ntypsd),ntypsy(natd)
      INTEGER, INTENT (IN) :: ulo_der(nlod)
      REAL,    INTENT (IN) :: aclo(nlod),bclo(nlod),cclo(nlod,nlod)
      REAL,    INTENT (IN) :: acnmt(0:lmaxd,nlod,nlhd)
      REAL,    INTENT (IN) :: bcnmt(0:lmaxd,nlod,nlhd)
      REAL,    INTENT (IN) :: ccnmt(nlod,nlod,nlhd)
      REAL,    INTENT (IN) :: dulon(nlod),uulon(nlod),vr(jmtd)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod),ello(nlod)
      REAL,    INTENT (IN) :: f(jmtd,2,0:lmaxd),g(jmtd,2,0:lmaxd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod)
      REAL,    INTENT (INOUT) :: qmtllo(0:lmaxd)
      REAL,    INTENT (INOUT) :: rho(jmtd,0:nlhd)
C     ..
C     .. Local Scalars ..
      REAL dsdum,usdum,c,c_1,c_2
      INTEGER j,l,lh,lo,lop,lp,nodedum
      REAL dus,ddn
C     ..
C     .. Local Arrays ..
      REAL,    ALLOCATABLE :: flo(:,:,:),glo(:,:)
      REAL filo(jmtd,2)
C     ..
      c = c_light(1.0)
      c_1 = 1.0 / neq
      c_2 = 1.0 /(neq*sfp)
c
      DO lo = 1,nlo
         l = llo(lo)
         qmtllo(l) = qmtllo(l) + (aclo(lo)*uulon(lo) +
     +                            bclo(lo)*dulon(lo)) * c_1
         DO lop = 1,nlo
            IF (llo(lop).EQ.l) THEN
               qmtllo(l) = qmtllo(l) + (cclo(lop,lo) * 
     +                              uloulopn(lop,lo)) * c_1
            END IF
         END DO
      END DO
      ALLOCATE ( flo(jmtd,2,nlod),glo(jmtd,2) )

c---> calculate the local ortital radial functions

      DO lo = 1,nlo
         l = llo(lo)
         CALL radsra(
     >               ello(lo),l,vr,r0,dx,jri,jmtd,c,
     <               usdum,dus,nodedum,flo(1,1,lo),flo(1,2,lo))
c+apw+lo
         IF (l_dulo(lo).or.ulo_der(lo).ge.1) THEN
c--->    calculate orthogonal energy derivative at e
           j = ulo_der(lo)
           IF(l_dulo(lo)) j = 1
           CALL radsrdn(
     >              ello(lo),l,vr,r0,dx,jri,jmtd,c,
     <              usdum,dsdum,ddn,nodedum,glo,filo, ! filo is a dummy array
     >              flo(1,1,lo),dus,j)
           DO j=1,jri
             flo(j,1,lo) = glo(j,1)
             flo(j,2,lo) = glo(j,2)
           ENDDO
           ddn = sqrt(ddn)
           IF(l_dulo(lo)) ddn=1.0
           flo(:,:,lo) = flo(:,:,lo)/ddn ! Normalize ulo (flo) if APW+lo is not used
         ENDIF
c-apw+lo
      END DO

c---> add the contribution of the local orbitals and flapw - lo cross-
c---> terms to the spherical chargedensity inside the muffin tins.

      DO lo = 1,nlo
         l = llo(lo)
         DO j = 1,jri
            rho(j,0) = rho(j,0) + c_2 *
     +                (aclo(lo) * ( f(j,1,l)*flo(j,1,lo) +
     +                              f(j,2,l)*flo(j,2,lo) ) +
     +                 bclo(lo) * ( g(j,1,l)*flo(j,1,lo) +
     +                              g(j,2,l)*flo(j,2,lo) ) )
         END DO
         DO lop = 1,nlo
            IF (llo(lop).EQ.l) THEN
               DO j = 1,jri
                  rho(j,0) = rho(j,0) + c_2 * cclo(lop,lo) *
     +                          ( flo(j,1,lop)*flo(j,1,lo) +
     +                            flo(j,2,lop)*flo(j,2,lo) )
               END DO
            END IF
         END DO
      END DO

c---> add the contribution of the local orbitals and flapw - lo cross-
c---> terms to the non-spherical chargedensity inside the muffin tins.

      DO lh = 1,nlh(ntypsy(nat))
         DO lp = 0,lmax
            DO lo = 1,nlo
               DO j = 1,jri
                  rho(j,lh) = rho(j,lh) + c_1 * (
     +                      acnmt(lp,lo,lh) * (f(j,1,lp)*flo(j,1,lo) +
     +                                         f(j,2,lp)*flo(j,2,lo) ) +
     +                      bcnmt(lp,lo,lh) * (g(j,1,lp)*flo(j,1,lo) +
     +                                         g(j,2,lp)*flo(j,2,lo) ) )
               END DO
            END DO
         END DO
         DO lo = 1,nlo
            DO lop = 1,nlo
               DO j = 1,jri
                  rho(j,lh) = rho(j,lh) + c_1 * ccnmt(lop,lo,lh) *
     +                                ( flo(j,1,lop)*flo(j,1,lo) +
     +                                  flo(j,2,lop)*flo(j,2,lo) )
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (flo,glo)

      END SUBROUTINE rhosphnlo
      END MODULE m_rhosphnlo
