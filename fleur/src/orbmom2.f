      MODULE m_orbmom2
c     ***************************************************************
c     perform the sum over m (for each l) and calculate the
c     spherical contribution to orbital moment.                
c     ***************************************************************
c
      CONTAINS
      SUBROUTINE orbmom2(
     >                   lmaxd,lmax,neq,itype,nlod,llod,nlo,llo,ddn,orb, 
     >                   uulon,dulon,uloulopn,orbl,orblo,
     <                   clmom)

      USE m_types, ONLY : t_orb,t_orbl,t_orblo
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,nlod,llod
      INTEGER, INTENT (IN) :: lmax,neq,itype,nlo
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod)
      REAL,    INTENT (IN) :: ddn(0:lmaxd),uulon(nlod),dulon(nlod)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod)
      TYPE (t_orb),  INTENT (IN) :: orb(0:lmaxd,-lmaxd:lmaxd)
      TYPE (t_orbl), INTENT (IN) :: orbl(nlod,-llod:llod)
      TYPE (t_orblo),INTENT (IN) :: orblo(nlod,nlod,-llod:llod)
      REAL,    INTENT (OUT) :: clmom(3)
C     ..
C     .. Local Scalars ..
      INTEGER l, m, ilo, ilop
      REAL qmtt, qmttx, qmtty, sumlm
      COMPLEX orbp, orbm
C     ..
C     .. Local Arrays ..
      REAL qmtl(0:lmaxd),qmtlx(0:lmaxd),qmtly(0:lmaxd)

      qmtt = 0.
      qmttx = 0.
      qmtty = 0.
      DO l = 0,lmax
c--->    lm-decomposed density for each atom type
         qmtl(l) = 0.
         qmtlx(l) = 0.
         qmtly(l) = 0.
         DO m = -l,l
c lz
            sumlm = m * (orb(l,m)%uu + orb(l,m)%dd * ddn(l) ) 
c lx,ly
            orbp = sqrt(real((l-m)*(l+m+1))) *
     *                ( orb(l,m)%uup + orb(l,m)%ddp * ddn(l) ) 

            orbm = sqrt(real((l+m)*(l-m+1))) *
     *                ( orb(l,m)%uum + orb(l,m)%ddm * ddn(l) )
c+gu
            IF (m.EQ.l)  orbp = cmplx(0.0,0.0)
            IF (m.EQ.-l) orbm = cmplx(0.0,0.0)
c+gu
            qmtl(l)  = qmtl(l)  + sumlm
            qmtlx(l) = qmtlx(l) + 0.5*( real(orbp)+ real(orbm))
            qmtly(l) = qmtly(l) + 0.5*(aimag(orbp)-aimag(orbm))
c 
         ENDDO
      ENDDO
!
! --> LO contribution
      DO ilo = 1, nlo
         l = llo(ilo)
         DO m = -l,l
            sumlm = m * (orbl(ilo,m)%uulo * uulon(ilo) +
     +                   orbl(ilo,m)%dulo * dulon(ilo) )

            orbp = sqrt(real((l-m)*(l+m+1))) *
     *                ( orbl(ilo,m)%uulop * uulon(ilo) +
     +                  orbl(ilo,m)%dulop * dulon(ilo) )

            orbm = sqrt(real((l+m)*(l-m+1))) *
     *                ( orbl(ilo,m)%uulom * uulon(ilo) +
     +                  orbl(ilo,m)%dulom * dulon(ilo) )

            IF (m.EQ.l)  orbp = cmplx(0.0,0.0)
            IF (m.EQ.-l) orbm = cmplx(0.0,0.0)

            qmtl(l)  = qmtl(l)  + sumlm
            qmtlx(l) = qmtlx(l) + 0.5*( real(orbp)+ real(orbm))
            qmtly(l) = qmtly(l) + 0.5*(aimag(orbp)-aimag(orbm))
         ENDDO
         DO ilop = 1, nlo
           IF (llo(ilop).EQ.l) THEN
             DO m = -l,l
               sumlm = m * orblo(ilo,ilop,m)%z * uloulopn(ilo,ilop)
               orbp = sqrt(real((l-m)*(l+m+1))) *
     *                     orblo(ilo,ilop,m)%p * uloulopn(ilo,ilop)
               orbm = sqrt(real((l+m)*(l-m+1))) *
     *                     orblo(ilo,ilop,m)%m * uloulopn(ilo,ilop)
               IF (m.EQ.l)  orbp = cmplx(0.0,0.0)
               IF (m.EQ.-l) orbm = cmplx(0.0,0.0)

               qmtl(l)  = qmtl(l)  + sumlm
               qmtlx(l) = qmtlx(l) + 0.5*( real(orbp)+ real(orbm))
               qmtly(l) = qmtly(l) + 0.5*(aimag(orbp)-aimag(orbm))
             ENDDO
           ENDIF
         ENDDO
      ENDDO
!
! --> sum up & print
      DO l = 0,lmax
         qmtl(l)  = qmtl(l)  / neq
         qmtlx(l) = qmtlx(l) / neq
         qmtly(l) = qmtly(l) / neq
         qmtt =  qmtt  + qmtl(l)
         qmttx = qmttx + qmtlx(l)
         qmtty = qmtty + qmtly(l)
      ENDDO
      clmom(1) = qmttx
      clmom(2) = qmtty
      clmom(3) = qmtt

      WRITE (6,FMT=8100) itype, (qmtl(l),l=0,3), qmtt
      WRITE (6,FMT=8100) itype, (qmtlx(l),l=0,3),qmttx
      WRITE (6,FMT=8100) itype, (qmtly(l),l=0,3),qmtty
 8100 FORMAT (' -->',i2,2x,4f9.5,2x,f9.5)

      END SUBROUTINE orbmom2
      END MODULE m_orbmom2
