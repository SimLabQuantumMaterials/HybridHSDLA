      MODULE m_forcea21lo
      CONTAINS
      SUBROUTINE force_a21_lo(
     >                        nobd,natd,ntypd,neigd,nlod,llod,lmd,
     >                        loplod,itype,natom,neq,lmax,we,eig,ne,
     >                        nlo,llo,nlol,lo1l,indmat,
     >                        acof,bcof,ccof,aveccof,bveccof,cveccof,
     >                        tuulo,tdulo,tuloulo,uulon,dulon,uloulopn,
     X                        a21)
c
c***********************************************************************
c This subroutine calculates the local orbital contribution to A21,
c which is the combination of the terms A17 and A20 according to the
c paper of R.Yu et al. (PRB vol.43 no.8 p.6411 1991).
c p.kurz nov. 1997
c***********************************************************************
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nobd,natd,ntypd,neigd,nlod,llod
      INTEGER, INTENT (IN) :: lmd,loplod,itype,natom,ne
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd)
      INTEGER, INTENT (IN) :: indmat(0:lmd,0:lmd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      INTEGER, INTENT (IN) :: nlol(0:llod,ntypd),lo1l(0:llod,ntypd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd),uulon(nlod,ntypd)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN) :: we(nobd),eig(neigd)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (IN) :: aveccof(3,nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bveccof(3,nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: cveccof(3,-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (IN) :: tuulo(0:lmd,-llod:llod,nlod)
      COMPLEX, INTENT (IN) :: tdulo(0:lmd,-llod:llod,nlod)
      COMPLEX, INTENT (IN) :: tuloulo(-llod:llod,-llod:llod,loplod)
      REAL, INTENT (INOUT) :: a21(3,natd)
C     ..
C     .. Local Scalars ..
      COMPLEX utulo,dtulo,cutulo,cdtulo,ulotulo
      INTEGER lo,lop,l,lp,m,mp,lm,lmp,iatom,ie,i,lolop,loplo,in
C     ..
C     ..
c*************** ABBREVIATIONS *****************************************
c ccof       : coefficient of the local orbital function (u_lo*Y_lm)
c cveccof    : is defined equivalently to aveccof, but with the LO-fct.
c tuulo,tdulo and tuloulo are the MT hamiltonian matrix elements of the
c local orbitals with the flapw basisfct. and with themselves.
c for information on nlo,llo,nlol,lo1l,uulon,dulon, and uloulopn see
c comments in setlomap.
c***********************************************************************

      DO lo = 1,nlo(itype)
         l = llo(lo,itype)
         DO m = -l,l
            lm = l* (l+1) + m
            DO lp = 0,lmax(itype)
               DO mp = -lp,lp
                  lmp = lp* (lp+1) + mp
                  DO iatom = natom,natom + neq(itype) - 1
c
c--->             check whether the t-matrixelement is 0
c--->             (indmat.EQ.-9999)
c
                     in = indmat(lmp,lm)
                     IF ((in.NE.-9999).OR.(lmp.EQ.lm)) THEN
                        utulo = tuulo(lmp,m,lo)
                        dtulo = tdulo(lmp,m,lo)
                        cutulo = conjg(tuulo(lmp,m,lo))
                        cdtulo = conjg(tdulo(lmp,m,lo))
                        DO ie = 1,ne
                           DO i = 1,3
                              a21(i,iatom)=a21(i,iatom)+2.0*aimag(
     +                             conjg(acof(ie,lmp,iatom))*utulo
     +                             *cveccof(i,m,ie,lo,iatom)
     +                           + conjg(bcof(ie,lmp,iatom))*dtulo
     +                             *cveccof(i,m,ie,lo,iatom)
     +                           + conjg(ccof(m,ie,lo,iatom))
     +                             *cutulo*aveccof(i,ie,lmp,iatom)
     +                           + conjg(ccof(m,ie,lo,iatom))
     +                             *cdtulo*bveccof(i,ie,lmp,iatom)
     +                                   )*we(ie)/neq(itype)
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO

               ENDDO
            ENDDO
            DO lop = 1,nlo(itype)
               lp = llo(lop,itype)
               DO mp = -lp,lp
                  lmp = lp* (lp+1) + mp
                  DO iatom = natom,natom + neq(itype) - 1
                     in = indmat(lmp,lm)
                     IF ((in.NE.-9999).OR.(lmp.EQ.lm)) THEN
                        IF (lo.GE.lop) THEN
                           lolop = (lo-1)*lo/2 + lop
                           ulotulo = tuloulo(m,mp,lolop)
                        ELSE
                           loplo = (lop-1)*lop/2 + lo
                           ulotulo = conjg(tuloulo(mp,m,loplo))
                        ENDIF
                        DO ie = 1,ne
                           DO i = 1,3
                              a21(i,iatom)=a21(i,iatom)+2.0*aimag(
     +                           + conjg(ccof(m,ie,lo,iatom))
     +                             *ulotulo*cveccof(i,mp,ie,lop,iatom)
     +                                   )*we(ie)/neq(itype)
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            DO iatom = natom,natom + neq(itype) - 1
               DO ie = 1,ne
                  DO i = 1,3
                     a21(i,iatom)=a21(i,iatom)
     -                           -2.0*aimag(
     +                            (
     +         conjg(acof(ie,lm,iatom))*cveccof(i,m,ie,lo,iatom)+
     +         conjg(ccof(m,ie,lo,iatom))*aveccof(i,ie,lm,iatom)
     +                            )*uulon(lo,itype)+

     +                            (
     +         conjg(bcof(ie,lm,iatom))*cveccof(i,m,ie,lo,iatom)+
     +         conjg(ccof(m,ie,lo,iatom))*bveccof(i,ie,lm,iatom)
     +                            )*dulon(lo,itype)
     +                            )*eig(ie)*we(ie)/neq(itype)
                  ENDDO
               ENDDO
            ENDDO
c--->       consider only the lop with l_lop = l_lo
            DO lop = lo1l(l,itype),(lo1l(l,itype)+nlol(l,itype)-1)
               DO iatom = natom,natom + neq(itype) - 1
                  DO ie = 1,ne
                     DO i = 1,3
                        a21(i,iatom)=a21(i,iatom)-2.0*aimag(
     +                       conjg(ccof(m,ie,lo,iatom))*
     +                       cveccof(i,m,ie,lop,iatom)*
     +                       uloulopn(lo,lop,itype)
     +                            )*eig(ie)*we(ie)/neq(itype)

                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
c--->    end of m loop
         ENDDO
c---> end of lo loop
      ENDDO

      END SUBROUTINE force_a21_lo
      END MODULE m_forcea21lo
