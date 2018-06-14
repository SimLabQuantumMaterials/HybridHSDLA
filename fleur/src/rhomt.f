      MODULE m_rhomt
      CONTAINS
      SUBROUTINE rhomt(
     >                 lmaxd,nobd,natd,ntypd,
     >                 we,ne,ntype,neq,lmax,acof,bcof,
     X                 uu,dd,du)
c     ***************************************************************
c     perform the sum over m (for each l) and bands to set up the
c     coefficient of spherical charge densities in subroutine
c     cdnval                                   c.l.fu
c     ***************************************************************

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,nobd,natd,ntypd
      INTEGER, INTENT (IN) :: ne,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmaxd* (lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmaxd* (lmaxd+2),natd)
      REAL,    INTENT (IN) :: we(nobd)
      REAL, INTENT (INOUT) :: dd(0:lmaxd,ntypd)
      REAL, INTENT (INOUT) :: du(0:lmaxd,ntypd)
      REAL, INTENT (INOUT) :: uu(0:lmaxd,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER i,l,lm,m,n,na,natom
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,real
C     ..
      natom = 0
      DO n = 1,ntype
         DO na = 1,neq(n)
            natom = natom + 1
            DO l = 0,lmax(n)
c     -----> sum over m
               DO m = -l,l
                  lm = l* (l+1) + m
c     -----> sum over occupied bands
                  DO i = 1,ne
                     uu(l,n) = uu(l,n) + we(i)*
     +                         real(acof(i,lm,natom)*conjg(acof(i,lm,
     +                         natom)))
                     dd(l,n) = dd(l,n) + we(i)*
     +                         real(bcof(i,lm,natom)*conjg(bcof(i,lm,
     +                         natom)))
                     du(l,n) = du(l,n) + we(i)*
     +                         real(acof(i,lm,natom)*conjg(bcof(i,lm,
     +                         natom)))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE rhomt
      END MODULE m_rhomt
