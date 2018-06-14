      MODULE m_int21
c-----------------------------------------------------------
c
c Integrates f(:,:,l,s) and g(:,:,l,s) (s=1,2) for given l
c and itype with different spins (s).
c Output is ..n21(l,itype), where .. is a (u,d) combination 
c dependet on the (f,g) combination used. 
c
c-----------------------------------------------------------
      CONTAINS
      SUBROUTINE int_21(
     >                  f,g,rmsh,dx,jri,jmtd,lmaxd,jspd,l,
     <                  uun21,udn21,dun21,ddn21)

      USE m_intgr, ONLY : intgr3
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,lmaxd,jmtd,jri,l
      REAL,    INTENT (IN) :: dx
      REAL,    INTENT (OUT):: uun21,udn21,dun21,ddn21
c     ... Array Arguments
      REAL,    INTENT (IN) :: rmsh(jri)
      REAL,    INTENT (IN) :: f(jmtd,2,0:lmaxd,jspd)
      REAL,    INTENT (IN) :: g(jmtd,2,0:lmaxd,jspd)
c     ...local scalars
      INTEGER iri
C     ...local arrays
      REAL        uu_tmp(jri)

      DO iri = 1, jri
         uu_tmp(iri) = f(iri,1,l,2)*f(iri,1,l,1)
     +               + f(iri,2,l,2)*f(iri,2,l,1)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            uun21)

      DO iri = 1, jri
         uu_tmp(iri) = f(iri,1,l,2)*g(iri,1,l,1)
     +               + f(iri,2,l,2)*g(iri,2,l,1)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            udn21)

      DO iri = 1, jri
         uu_tmp(iri) = g(iri,1,l,2)*f(iri,1,l,1)
     +               + g(iri,2,l,2)*f(iri,2,l,1)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            dun21)

      DO iri = 1, jri
         uu_tmp(iri) = g(iri,1,l,2)*g(iri,1,l,1)
     +               + g(iri,2,l,2)*g(iri,2,l,1)
      ENDDO
      CALL intgr3(uu_tmp,rmsh,dx,jri,
     <            ddn21)

      END SUBROUTINE int_21
      END MODULE m_int21
