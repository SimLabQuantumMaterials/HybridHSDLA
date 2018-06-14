      MODULE m_forcea21U
      CONTAINS
      SUBROUTINE force_a21_U(
     >                       nobd,natd,nlod,llod,lmd,lmaxb,
     >                       natom,neq,lmax,we,ne,
     >                       nlo,llo,lda_u,ddn,v_mmp,
     >                       acof,bcof,ccof,aveccof,bveccof,cveccof,
     >                       uulon,dulon,
     X                       a21)
!
!***********************************************************************
! This subroutine calculates the lda+U contribution to the HF forces, 
! similar to the A21 term, according to eqn. (22) of F. Tran et al.
! Comp.Phys.Comm. 179 (2008) 784-790
!***********************************************************************
!
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nobd,natd,nlod,llod
      INTEGER, INTENT (IN) :: lmd,natom,ne,lmaxb
      INTEGER, INTENT (IN) :: lmax,lda_u,neq,nlo
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod)
      REAL,    INTENT (IN) :: dulon(nlod),uulon(nlod)
      REAL,    INTENT (IN) :: we(nobd),ddn(0:lmax)
      COMPLEX, INTENT (IN) :: v_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (IN) :: aveccof(3,nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bveccof(3,nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: cveccof(3,-llod:llod,nobd,nlod,natd)
      REAL, INTENT (INOUT) :: a21(3,natd)
C     ..
C     .. Local Scalars ..
      COMPLEX v_a,v_b,v_c,p1,p2,p3
      INTEGER lo,lop,l,lp,m,mp,lm,lmp,iatom,ie,i
C     ..
C     ..
c*************** ABBREVIATIONS *****************************************
c ccof       : coefficient of the local orbital function (u_lo*Y_lm)
c cveccof    : is defined equivalently to aveccof, but with the LO-fct.
c for information on nlo,llo,uulon,dulon, and uloulopn see
c comments in setlomap.
c***********************************************************************

      IF (lda_u.GE.0) THEN
        l = lda_u
!
! Add contribution for the regular LAPWs (like force_a21, but with
! the potential matrix, v_mmp, instead of the tuu, tdd ...)
!
        DO m = -l,l
          lm = l* (l+1) + m
          DO mp = -l,l
            lmp = l* (l+1) + mp
            v_a = v_mmp(m,mp) 
            v_b = v_mmp(m,mp) * ddn(l) 

            DO iatom = natom,natom + neq - 1
              DO ie = 1,ne
                DO i = 1,3

                  p1 = ( conjg(acof(ie,lm,iatom)) * v_b )
     +                    * aveccof(i,ie,lmp,iatom)
                  p2 = ( conjg(bcof(ie,lm,iatom)) * v_b )
     +                    * bveccof(i,ie,lmp,iatom) 
                  a21(i,iatom) = a21(i,iatom) + 2.0*aimag(
     +                        p1 + p2 ) *we(ie)/neq

! no idea, why this did not work with ifort:
!                  a21(i,iatom) = a21(i,iatom) + 2.0*aimag(
!     +                         conjg(acof(ie,lm,iatom)) * v_a *
!     +                         *aveccof(i,ie,lmp,iatom)   +
!     +                         conjg(bcof(ie,lm,iatom)) * v_b *
!     +                         *bveccof(i,ie,lmp,iatom)   )
!     +                                       *we(ie)/neq
                ENDDO
              ENDDO
            ENDDO

          ENDDO ! mp
        ENDDO   ! m
!
! If there are also LOs on this atom, with the same l as
! the one of LDA+U, add another few terms
!
        DO lo = 1,nlo
          l = llo(lo)
          IF ( l == lda_u ) THEN

            DO m = -l,l
              lm = l* (l+1) + m
              DO mp = -l,l
                lmp = l* (l+1) + mp
                v_a = v_mmp(m,mp)
                v_b = v_mmp(m,mp) * uulon(lo)
                v_c = v_mmp(m,mp) * dulon(lo)

                DO iatom = natom,natom + neq - 1
                  DO ie = 1,ne
                    DO i = 1,3

                    p1 = v_a * ( conjg(ccof(m,ie,lo,iatom)) 
     +                         * cveccof(i,mp,ie,lo,iatom) )
                    p2 = v_b * ( conjg(acof(ie,lm,iatom))
     +                       * cveccof(i,mp,ie,lo,iatom) +
     +                         conjg(ccof(m,ie,lo,iatom))
     +                       *   aveccof(i,ie,lmp,iatom) )
                    p3 = v_c * ( conjg(bcof(ie,lm,iatom))
     +                       * cveccof(i,mp,ie,lo,iatom) +
     +                         conjg(ccof(m,ie,lo,iatom))
     +                       *   bveccof(i,ie,lmp,iatom) )
                    a21(i,iatom) = a21(i,iatom) + 2.0*aimag(
     +                      p1 + p2 + p3 )*we(ie)/neq

                    ENDDO
                  ENDDO
                ENDDO

              ENDDO
            ENDDO

          ENDIF   ! l == lda_u
        ENDDO     ! lo = 1,nlo

      ENDIF

      END SUBROUTINE force_a21_U
      END MODULE m_forcea21U
