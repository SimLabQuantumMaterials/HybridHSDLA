      MODULE m_rhomtlo
c
c***********************************************************************
c This subroutine is the equivalent of rhomt for the local orbital
c contributions to the charge.
c aclo,bclo,cclo are the equivalents of uu,ud,dd in rhomt
c p.kurz sept. 1996
c***********************************************************************
c
      CONTAINS
      SUBROUTINE rhomtlo(
     >                   nobd,natd,lmd,nlod,llod,ntypd,
     >                   ntype,neq,ne,we,acof,bcof,ccof,nlo,llo,invsat,
     X                   aclo,bclo,cclo)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: natd,lmd,nlod,llod,nobd,ntypd
      INTEGER, INTENT (IN) :: ne,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),invsat(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: we(nobd)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd)
      REAL,    INTENT (INOUT):: aclo(nlod,ntypd),bclo(nlod,ntypd)
      REAL,    INTENT (INOUT):: cclo(nlod,nlod,ntypd)
C     ..
C     .. Local Scalars ..
      REAL invsfct
      INTEGER i,l,lm,lo,lop,m,natom,nn,ntyp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,real
C     ..

      natom = 0
c---> loop over atoms
      DO ntyp = 1,ntype
         DO nn = 1,neq(ntyp)
            natom = natom + 1
            IF ((invsat(natom).EQ.0) .OR. (invsat(natom).EQ.1)) THEN
               IF (invsat(natom).EQ.0) invsfct = 1.0
               IF (invsat(natom).EQ.1) invsfct = 2.0
c--->       loop over the local orbitals
               DO lo = 1,nlo(ntyp)
                  l = llo(lo,ntyp)
c--->          contribution of cross terms flapw - local orbitals
                  DO m = -l,l
                     lm = l* (l+1) + m
                     DO i = 1,ne
                        aclo(lo,ntyp) = aclo(lo,ntyp) +
     +                                  we(i)*invsfct*2*real(conjg(acof
     +                                 (i,lm,natom))*ccof(m,i,lo,natom))
                        bclo(lo,ntyp) = bclo(lo,ntyp) +
     +                                  we(i)*invsfct*2*real(conjg(bcof
     +                                 (i,lm,natom))*ccof(m,i,lo,natom))
                     END DO
                  END DO
c--->          contribution of local orbital - local orbital terms
c--->          loop over lo'
                  DO lop = 1,nlo(ntyp)
                     IF (llo(lop,ntyp).EQ.l) THEN
                        DO m = -l,l
                           DO i = 1,ne
                              cclo(lop,lo,ntyp) = cclo(lop,lo,ntyp) +
     +                                            we(i)*invsfct*
     +                                            real(conjg(ccof(m,i,
     +                                            lop,natom))*
     +                                            ccof(m,i,lo,natom))
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO


      END SUBROUTINE rhomtlo
      END MODULE m_rhomtlo
