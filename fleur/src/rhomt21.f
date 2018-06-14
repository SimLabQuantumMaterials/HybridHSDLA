      MODULE m_rhomt21
c     ***************************************************************
c     perform the sum over m (for each l) and bands to set up the
c     coefficient of spherical charge densities in subroutine 
c     cdnval                                 
c     for offdiagonal matrix-elements in case of noncollinear magnetism 
c     FF
c     ***************************************************************
      CONTAINS
      SUBROUTINE rhomt21(
     >                   lmaxd,nobd,natd,ntypd,jspd,
     >                   we,ne,ntype,neq,lmax,acof,bcof,
     >                   nlod,llod,nlo,llo,ccof,
     X                   mt21,lo21,uloulop21)

      USE m_types, ONLY : t_mt21,t_lo21
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,nobd,natd,ntypd,jspd,nlod,llod
      INTEGER, INTENT (IN) :: ne,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmaxd* (lmaxd+2),natd,jspd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmaxd* (lmaxd+2),natd,jspd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd,jspd)
      REAL,    INTENT (IN) :: we(nobd)
      TYPE (t_mt21), INTENT (INOUT) :: mt21(0:lmaxd,ntypd)
      TYPE (t_lo21), INTENT (INOUT) :: lo21(nlod,ntypd)
      COMPLEX,       INTENT (INOUT) :: uloulop21(nlod,nlod,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER   i,l,lm,m,itype,na,natom,lo,lop

C     .. Intrinsic Functions ..
      INTRINSIC conjg,real
C     ..
      natom = 0
      DO itype = 1,ntype
         DO na = 1,neq(itype)
            natom = natom + 1
!
!--->       normal u, du contribution
!
            DO l = 0,lmax(itype)
              DO m = -l,l
                lm = l* (l+1) + m
c--->           sum over occupied bands
                DO i = 1,ne
                  mt21(l,itype)%uu = mt21(l,itype)%uu + we(i)*
     +                      conjg(acof(i,lm,natom,2))*acof(i,lm,natom,1)
                  mt21(l,itype)%ud = mt21(l,itype)%ud + we(i)*
     +                      conjg(acof(i,lm,natom,2))*bcof(i,lm,natom,1)
                  mt21(l,itype)%du = mt21(l,itype)%du + we(i)*
     +                      conjg(bcof(i,lm,natom,2))*acof(i,lm,natom,1)
                  mt21(l,itype)%dd = mt21(l,itype)%dd + we(i)*
     +                      conjg(bcof(i,lm,natom,2))*bcof(i,lm,natom,1)
                ENDDO ! i = 1,ne
              ENDDO   ! m = -l,l
            ENDDO     ! l
!
!--->       loop over the local orbitals
!
            DO lo = 1,nlo(itype)
              l = llo(lo,itype)
c--->         contribution of cross terms flapw - local orbitals
              DO m = -l,l
                lm = l* (l+1) + m
                DO i = 1,ne
                  lo21(lo,itype)%uulo = lo21(lo,itype)%uulo + we(i)*
     +                    conjg(acof(i,lm,natom,2))*ccof(m,i,lo,natom,1)
                  lo21(lo,itype)%dulo = lo21(lo,itype)%dulo + we(i)*
     +                    conjg(bcof(i,lm,natom,2))*ccof(m,i,lo,natom,1)
                  lo21(lo,itype)%ulou = lo21(lo,itype)%ulou + we(i)*
     +                    conjg(acof(i,lm,natom,1))*ccof(m,i,lo,natom,2)
                  lo21(lo,itype)%ulod = lo21(lo,itype)%ulod + we(i)*
     +                    conjg(bcof(i,lm,natom,1))*ccof(m,i,lo,natom,2)
                ENDDO
              ENDDO
c--->         contribution of local orbital - local orbital terms
c--->         loop over lo'
              DO lop = 1,nlo(itype)
                IF (llo(lop,itype).EQ.l) THEN
                  DO m = -l,l
                    DO i = 1,ne
                      uloulop21(lop,lo,itype) = uloulop21(lop,lo,itype)+
     +                               we(i)*conjg(ccof(m,i,lop,natom,2))*
     +                                           ccof(m,i,lo, natom,1)
                    ENDDO ! i = 1,ne
                  ENDDO   ! m = -l,l
                ENDIF
              ENDDO     ! lop
            ENDDO       ! lo

         ENDDO          ! na
      ENDDO             ! itype

      END SUBROUTINE rhomt21
      END MODULE m_rhomt21
