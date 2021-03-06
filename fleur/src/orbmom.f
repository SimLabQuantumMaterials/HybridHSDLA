      MODULE m_orbmom
!     ***************************************************************
!     perform the sum over m (for each l) and bands to set up the
!     coefficient of spherical contribution to orbital moment.
!     all quantities are in the local spin-frame
!     ***************************************************************

      CONTAINS
      SUBROUTINE orbmom(
     >                  ntypd,lmaxd,natd,nobd,lmd,
     >                  ntype,neq,ne,lmax,we,acof,bcof,
     >                  nlod,llod,nlo,llo,ccof,
     X                  orb,orbl,orblo)

      USE m_types, ONLY : t_orb,t_orbl,t_orblo
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,lmaxd,natd,nobd,lmd,nlod,llod
      INTEGER, INTENT (IN) :: ne,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lmax(ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd)
      REAL,    INTENT (IN) :: we(nobd)
      TYPE (t_orb),  INTENT (INOUT) :: orb(0:lmaxd,-lmaxd:lmaxd,ntypd)
      TYPE (t_orbl), INTENT (INOUT) :: orbl(nlod,-llod:llod,ntypd)
      TYPE (t_orblo),INTENT (INOUT) :: orblo(nlod,nlod,-llod:llod,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER i,l,lm,m,n,na,natom,ilo,ilop
      COMPLEX czero
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,cmplx
C     ..
      czero = cmplx(0.0,0.0)

      natom = 0
      DO n = 1,ntype
         DO na = 1,neq(n)
            natom = natom + 1

            DO 30 l = 0,lmax(n)
c     -----> sum over m
               DO 20 m = -l,l
                  lm = l* (l+1) + m
c     -----> sum over occupied bands
                  DO 10 i = 1,ne
c coeff. for lz ->
                      orb(l,m,n)%uu = orb(l,m,n)%uu +
     +                               we(i)*acof(i,lm,natom)*
     +                               conjg(acof(i,lm,natom))
                      orb(l,m,n)%dd = orb(l,m,n)%dd +
     +                               we(i)*bcof(i,lm,natom)*
     +                               conjg(bcof(i,lm,natom))
c coeff. for l+ <M'|l+|M> with respect to M ->
                      IF (m.NE.l) THEN
                          orb(l,m,n)%uup = orb(l,m,n)%uup +
     +                                   we(i)*acof(i,lm,natom)*
     +                                   conjg(acof(i,lm+1,natom))
                          orb(l,m,n)%ddp = orb(l,m,n)%ddp +
     +                                   we(i)*bcof(i,lm,natom)*
     +                                   conjg(bcof(i,lm+1,natom))
                       ELSE
                          orb(l,m,n)%uup = czero
                          orb(l,m,n)%ddp = czero
                       ENDIF
c coeff. for l- <M'|l-|M> with respect to M ->
                       IF (m.NE.-l) THEN
                          orb(l,m,n)%uum = orb(l,m,n)%uum +
     +                                    we(i)*acof(i,lm,natom)*
     +                                    conjg(acof(i,lm-1,natom))
                          orb(l,m,n)%ddm = orb(l,m,n)%ddm +
     +                                    we(i)*bcof(i,lm,natom)*
     +                                    conjg(bcof(i,lm-1,natom))
                       ELSE
                          orb(l,m,n)%uum = czero
                          orb(l,m,n)%ddm = czero
                       ENDIF
   10             ENDDO
   20          ENDDO
   30       ENDDO
!
! --> Local Orbital contribution: u,lo part
!
            DO ilo = 1, nlo(n)
              l = llo(ilo,n)
              DO m = -l, l
                lm = l* (l+1) + m
                DO i = 1,ne
                  orbl(ilo,m,n)%uulo = orbl(ilo,m,n)%uulo + we(i) * (
     +                    acof(i,lm,natom)* conjg(ccof(m,i,ilo,natom)) +
     +                    ccof(m,i,ilo,natom)* conjg(acof(i,lm,natom)) )
                  orbl(ilo,m,n)%dulo = orbl(ilo,m,n)%dulo + we(i) * (
     +                    bcof(i,lm,natom)* conjg(ccof(m,i,ilo,natom)) +
     +                    ccof(m,i,ilo,natom)* conjg(bcof(i,lm,natom)) )
                  IF (m.NE.l) THEN
                    orbl(ilo,m,n)%uulop = orbl(ilo,m,n)%uulop + we(i) *(
     +                   acof(i,lm,natom)* conjg(ccof(m+1,i,ilo,natom))+
     +                   ccof(m,i,ilo,natom)* conjg(acof(i,lm+1,natom)))
                    orbl(ilo,m,n)%dulop = orbl(ilo,m,n)%dulop + we(i) *(
     +                   bcof(i,lm,natom)* conjg(ccof(m+1,i,ilo,natom))+
     +                   ccof(m,i,ilo,natom)* conjg(bcof(i,lm+1,natom)))
                  ELSE
                    orbl(ilo,m,n)%uulop = czero
                    orbl(ilo,m,n)%dulop = czero
                  ENDIF
                  IF (m.NE.-l) THEN
                    orbl(ilo,m,n)%uulom = orbl(ilo,m,n)%uulom + we(i) *(
     +                   acof(i,lm,natom)* conjg(ccof(m-1,i,ilo,natom))+
     +                   ccof(m,i,ilo,natom)* conjg(acof(i,lm-1,natom)))
                    orbl(ilo,m,n)%dulom = orbl(ilo,m,n)%dulom + we(i) *(
     +                   bcof(i,lm,natom)* conjg(ccof(m-1,i,ilo,natom))+
     +                   ccof(m,i,ilo,natom)* conjg(bcof(i,lm-1,natom)))
                  ELSE
                    orbl(ilo,m,n)%uulom = czero
                    orbl(ilo,m,n)%dulom = czero
                  ENDIF
                ENDDO  ! sum over eigenstates (i)
              ENDDO    ! loop over m
!
! --> lo,lo' part           
!
              DO ilop = 1, nlo(n)
                IF (llo(ilop,n).EQ.l) THEN
                  DO m = -l, l
                    DO i = 1,ne
                      orblo(ilo,ilop,m,n)%z = orblo(ilo,ilop,m,n)%z +
     +                               we(i) *   ccof(m,i,ilo, natom) *
     +                                  conjg( ccof(m,i,ilop,natom) ) 
                      IF (m.NE.l) THEN
                        orblo(ilo,ilop,m,n)%p = orblo(ilo,ilop,m,n)%p +
     +                                we(i) *  ccof(m,  i,ilo, natom) *
     +                                  conjg( ccof(m+1,i,ilop,natom) ) 
                      ELSE
                        orblo(ilo,ilop,m,n)%p = czero
                      ENDIF
                      IF (m.NE.-l) THEN
                        orblo(ilo,ilop,m,n)%m = orblo(ilo,ilop,m,n)%m +
     +                                we(i) *  ccof(m,  i,ilo, natom) *
     +                                  conjg( ccof(m-1,i,ilop,natom) )  
                      ELSE
                        orblo(ilo,ilop,m,n)%m = czero
                      ENDIF
                    ENDDO  ! sum over eigenstates (i)
                  ENDDO    ! loop over m
                ENDIF
              ENDDO      ! loop over lo's (ilop)

            ENDDO      ! loop over lo's (ilo)

         ENDDO ! sum over equiv atoms (na)
      ENDDO    ! loop over atom types (n)

      RETURN
      END SUBROUTINE orbmom
      END MODULE m_orbmom
