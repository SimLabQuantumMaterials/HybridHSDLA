      MODULE m_int21lo
!-----------------------------------------------------------
!
! Integrates combinations of flo(:,:,llo,s) with itself, f(:,:,l,s) 
! and g(:,:,l,s) (s=1,2) for given llo and itype with different 
! spins (s).
! Output is ..n21(l,itype), where .. is a combination of (u,d) and
! ulo dependent on the (f,g) combination used. Also ..n12 and 
! uloulopn21 are calculated.
!
!-----------------------------------------------------------
      CONTAINS
      SUBROUTINE int_21lo(
     >                    f,g,rmsh,dx,jri,jmtd,lmaxd,jspd,
     >                    flo,llo,ilo,nlo,nlod,
     <                    uulon21,dulon21,uulon12,dulon12,uloulopn21)

      USE m_intgr, ONLY : intgr3
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,jmtd,jspd,nlod,jri,nlo,ilo
      REAL,    INTENT (IN) :: dx
c     ... Array Arguments
      INTEGER, INTENT (IN) :: llo(nlo)
      REAL,    INTENT (IN) :: rmsh(jri)
      REAL,    INTENT (IN) :: f(jmtd,2,0:lmaxd,jspd)
      REAL,    INTENT (IN) :: g(jmtd,2,0:lmaxd,jspd)
      REAL,    INTENT (IN) :: flo(jmtd,2,nlod,jspd)
      REAL,    INTENT (OUT):: uulon21,uulon12
      REAL,    INTENT (OUT):: dulon21,dulon12
      REAL,    INTENT (OUT):: uloulopn21(nlod,nlod)

c     ...local scalars
      INTEGER iri,l,lp,ilop
C     ...local arrays
      REAL    uu_tmp(jri)

!
! --> norm of product of u and ulo:
!
      l = llo(ilo)
      DO iri = 1, jri
         uu_tmp(iri) = f(iri,1,l,2)*flo(iri,1,ilo,1)
     +               + f(iri,2,l,2)*flo(iri,2,ilo,1)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            uulon21)
      DO iri = 1, jri
         uu_tmp(iri) = f(iri,1,l,1)*flo(iri,1,ilo,2)
     +               + f(iri,2,l,1)*flo(iri,2,ilo,2)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            uulon12)
!
! --> norm of product of du and ulo:
!
      DO iri = 1, jri
         uu_tmp(iri) = g(iri,1,l,2)*flo(iri,1,ilo,1)
     +               + g(iri,2,l,2)*flo(iri,2,ilo,1)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            dulon21)
      DO iri = 1, jri
         uu_tmp(iri) = g(iri,1,l,1)*flo(iri,1,ilo,2)
     +               + g(iri,2,l,1)*flo(iri,2,ilo,2)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            dulon12)
!
! --> norm of product of ulo and ulo':
!
      DO ilop = 1, nlo
        lp = llo(ilop)
        IF (l.EQ.lp) THEN
          DO iri = 1, jri
            uu_tmp(iri) = flo(iri,1,ilo,2)*flo(iri,1,ilop,1)
     +                  + flo(iri,2,ilo,2)*flo(iri,2,ilop,1)
          ENDDO
          CALL intgr3(uu_tmp,rmsh,dx,jri,
     <                uloulopn21(ilo,ilop))
        ELSE
          uloulopn21(ilo,ilop) = 0.0
        ENDIF
      ENDDO

      END SUBROUTINE int_21lo
      END MODULE m_int21lo
