      MODULE m_rhonmtlo
c
c***********************************************************************
c This subroutine is the equivalent of rhomt for the local orbital
c contributions to the charge.
c acnmt, bcnmt, ccnmt are the equivalents of uunmt, ddnmt, udnmt dunmt
c in rhonmt
c p.kurz sept. 1996
c***********************************************************************
c
      CONTAINS
      SUBROUTINE rhonmtlo(
     >                   natd,lmd,nlod,memd,nobd,ntypd,nlhd,lmaxd,
     >                   ntype,neq,lmax,ne,we,acof,bcof,ccof,llod,
     >                   clnu,nmem,mlh,nlh,llh,ntypsy,ntypsd,nlo,llo,
     X                   acnmt,bcnmt,ccnmt)
      USE m_gaunt,ONLY:gaunt1
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: natd,lmd,nlod,llod,nobd,ntypd,nlhd,lmaxd
      INTEGER, INTENT (IN) :: ne,ntype,memd,ntypsd
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),ntypsy(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: we(nobd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd)
      REAL,    INTENT (INOUT) :: acnmt(0:lmaxd,nlod,nlhd,ntypd)
      REAL,    INTENT (INOUT) :: bcnmt(0:lmaxd,nlod,nlhd,ntypd)
      REAL,    INTENT (INOUT) :: ccnmt(nlod,nlod,nlhd,ntypd)
C     ..
C     .. Local Scalars ..
      COMPLEX ci,cmv,fact
      INTEGER i,jmem,l,lh,lmp,lo,lop,lp,lpmax,lpmax0,lpmin,lpmin0,
     +        lpp,m,mp,mpp,na,neqat0,nn,ntyp
C     ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,conjg,max,min,mod,real
C     ..

      ci = cmplx(0.0,1.0)

c---> for optimal performance consider only
c---> those combinations of l,l',l'',m,m',m'' that satisfy the three
c---> conditions for non-zero gaunt-coeff. i.e.
c---> |l - l''| <= l' <= l + l'' (triangular condition)
c---> m' + m'' = m and l + l' + l'' even

      neqat0 = 0
      DO ntyp = 1,ntype
c--->    loop over the lattice harmonics
         DO lh = 1,nlh(ntypsy(neqat0+1))
            lpp = llh(lh,ntypsy(neqat0+1))
            DO jmem = 1,nmem(lh,ntypsy(neqat0+1))
               mpp = mlh(jmem,lh,ntypsy(neqat0+1))
               cmv = conjg(clnu(jmem,lh,ntypsy(neqat0+1)))
               DO lo = 1,nlo(ntyp)
                  l = llo(lo,ntyp)
                  lpmin0 = abs(l-lpp)
                  lpmax0 = l + lpp
c--->             check that lpmax is smaller than the max l of the
c--->             wavefunction expansion at this atom
                  lpmax = min(lpmax0,lmax(ntyp))
c--->             make sure that l + l'' + lpmax is even
                  lpmax = lpmax - mod(l+lpp+lpmax,2)
                  DO m = -l,l

c--->                add flapw - local orbital cross-terms

c--->                add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
c--->                note that gaunt1(l,lp,lpp,m,mp,mpp) computes the
c--->                integral of conjg(y(l,m))*y(lp,mp)*y(lpp,mpp),
c--->                however, since the gaunt coef. are real, this is
c--->                the same as int. y(l,m)*conjg(y(lp,mp)*y(lpp,mpp))
                     mp = m - mpp
                     lpmin = max(lpmin0,abs(mp))
c--->                make sure that l + l'' + lpmin is even
                     lpmin = lpmin + mod(abs(lpmax-lpmin),2)
c--->                loop over l'
                     DO lp = lpmin,lpmax,2
                        lmp = lp* (lp+1) + mp
                        fact = cmv* (ci** (l-lp))*
     +                         gaunt1(l,lp,lpp,m,mp,mpp,lmaxd)
                        na = neqat0
                        DO nn = 1,neq(ntyp)
                           na = na + 1
                           DO i = 1,ne
                              acnmt(lp,lo,lh,ntyp) =
     +                        acnmt(lp,lo,lh,ntyp) + we(i) * real(  
     +                                fact * conjg(acof(i,lmp,na)) 
     +                                     *       ccof(m,i,lo,na))
                              bcnmt(lp,lo,lh,ntyp) =
     +                        bcnmt(lp,lo,lh,ntyp) + we(i) * real(  
     +                                fact * conjg(bcof(i,lmp,na)) 
     +                                     *       ccof(m,i,lo,na))
                           END DO
                        END DO
                     END DO

c--->                add terms containing gaunt1(lp,l,lpp,mp,m,mpp)
                     mp = m + mpp
                     lpmin = max(lpmin0,abs(mp))
c--->                make sure that l + l'' + lpmin is even
                     lpmin = lpmin + mod(abs(lpmax-lpmin),2)
c--->                loop over l'
                     DO lp = lpmin,lpmax,2
                        lmp = lp* (lp+1) + mp
                        fact = cmv* (ci** (lp-l))*
     +                         gaunt1(lp,l,lpp,mp,m,mpp,lmaxd)
                        na = neqat0
                        DO nn = 1,neq(ntyp)
                           na = na + 1
                           DO i = 1,ne
                              acnmt(lp,lo,lh,ntyp) =
     +                        acnmt(lp,lo,lh,ntyp) + we(i) * real(  
     +                                fact * conjg(ccof(m,i,lo,na)) 
     +                                     *       acof(i,lmp,na) )
!     +                                fact *      (ccof(m,i,lo,na)) 
!     +                                     * conjg(acof(i,lmp,na)))
                              bcnmt(lp,lo,lh,ntyp) =
     +                        bcnmt(lp,lo,lh,ntyp) + we(i) * real(  
     +                                fact * conjg(ccof(m,i,lo,na)) 
     +                                     *       bcof(i,lmp,na) )
!     +                                fact *      (ccof(m,i,lo,na)) 
!     +                                     * conjg(bcof(i,lmp,na)))
                           END DO
                        END DO
                     END DO

c--->                add local orbital - local orbital terms
                     DO lop = 1,nlo(ntyp)
                        lp = llo(lop,ntyp)

c--->                   add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
                        mp = m - mpp
                        IF ((abs(l-lpp).LE.lp) .AND.
     +                      (lp.LE. (l+lpp)) .AND.
     +                      (mod(l+lp+lpp,2).EQ.0) .AND.
     +                      (abs(mp).LE.lp)) THEN
                           fact = cmv* (ci** (l-lp))*
     +                            gaunt1(l,lp,lpp,m,mp,mpp,lmaxd)
                           na = neqat0
                           DO nn = 1,neq(ntyp)
                              na = na + 1
                              DO i = 1,ne
                                 ccnmt(lop,lo,lh,ntyp) =
     +                           ccnmt(lop,lo,lh,ntyp) + we(i) * real(
     +                                   fact * conjg(ccof(mp,i,lop,na))
     +                                        *       ccof(m ,i,lo ,na))
                              END DO
                           END DO
                        END IF

                     END DO
                  END DO
               END DO
            END DO
         END DO
         neqat0 = neqat0 + neq(ntyp)
      END DO

      END SUBROUTINE rhonmtlo
      END MODULE m_rhonmtlo
