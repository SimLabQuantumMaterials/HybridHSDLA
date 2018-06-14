      MODULE m_rhonmt
      CONTAINS
      SUBROUTINE rhonmt(
     >                  lmaxd,ntypd,memd,nlhd,ntypsd,natd,nobd,
     >                  we,ne,ntype,nsymt,neq,lmax,
     >                  clnu,nmem,mlh,nlh,llh,ntypsy,
     >                  acof,bcof,
     X                  uunmt,ddnmt,udnmt,dunmt)
c     *************************************************************
c     subroutine sets up the coefficients of non-sphereical
c     muffin-tin density                          c.l.fu
c     *************************************************************
      USE m_gaunt,ONLY:gaunt1
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT(IN) :: lmaxd,ntypd,memd,nlhd,ntypsd,natd,nobd
      INTEGER,INTENT(IN) :: ne,ntype,nsymt
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT(IN) :: neq(ntypd),ntypsy(natd),lmax(ntypd)
      INTEGER, INTENT(IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT(IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      COMPLEX, INTENT(IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT(IN) :: acof(nobd,0:lmaxd* (lmaxd+2),natd)
      COMPLEX, INTENT(IN) :: bcof(nobd,0:lmaxd* (lmaxd+2),natd)
      REAL,    INTENT(IN) :: we(nobd)
      REAL,INTENT(INOUT) :: ddnmt(0:(lmaxd* (lmaxd+3))/2,nlhd,ntypd)
      REAL,INTENT(INOUT) :: dunmt(0:(lmaxd* (lmaxd+3))/2,nlhd,ntypd)
      REAL,INTENT(INOUT) :: udnmt(0:(lmaxd* (lmaxd+3))/2,nlhd,ntypd)
      REAL,INTENT(INOUT) :: uunmt(0:(lmaxd* (lmaxd+3))/2,nlhd,ntypd)
C     ..
C     .. Local Scalars ..
      COMPLEX cconst,cil,cmv,ci
      REAL coef
      INTEGER jmem,l,lcond,lh,llp,llpmax,lm,lmp,lp,lphi,lplow,lplow0,lv,
     +        m,mp,mv,na,natom,nb,nn,ns,nt
C     ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,iabs,max,mod,real,cmplx
C     ..
      ci = cmplx(0.0,1.0)
c
      DO 100 ns = 1,nsymt
         DO 90 lh = 1,nlh(ns)
            lv = llh(lh,ns)
            DO 80 jmem = 1,nmem(lh,ns)
               mv = mlh(jmem,lh,ns)
               cmv = conjg(clnu(jmem,lh,ns))
               DO 70 l = 0,lmaxd
                  m_loop: DO m = -l,l
                     lm = l* (l+1) + m
                     mp = m - mv
c     -----> set up the lower and upper limit of lp in such a way that
c     -----> lp+l+lv is even, lp<=l, and (lp,l,lv) satisfies the
c     -----> triangular relation
                     lplow0 = iabs(l-lv)
                     lphi = l - mod(lv,2)
                     lplow = max(lplow0,iabs(mp))
                     lcond = iabs(lphi-lplow)
                     lplow = lplow + mod(lcond,2)
                     IF (lplow.GT.lphi) CYCLE m_loop
                     DO 50 lp = lplow,lphi,2
                        cil = ci** (l-lp)
                        lmp = lp* (lp+1) + mp
                        IF (lmp.GT.lm) CYCLE m_loop
                        llp = (l* (l+1))/2 + lp
c     -----> gaunt's coefficient
                        coef = 2.*gaunt1(l,lv,lp,m,mv,mp,lmaxd)
                        IF (lmp.EQ.lm) coef = coef/2.
                        cconst = coef* (cil*cmv)
                        natom = 0
                        DO 40 nn = 1,ntype
                           llpmax = (lmax(nn)* (lmax(nn)+3))/2
                           IF (llp.GT.llpmax) GO TO 30
                           nt = natom
                           DO 20 na = 1,neq(nn)
                              nt = nt + 1
                              IF (ntypsy(nt).EQ.ns) THEN
                                DO 10 nb = 1,ne
                                 uunmt(llp,lh,nn) = uunmt(llp,lh,nn) +
     +                                              we(nb)*real(cconst*
     +                                              acof(nb,lm,nt)*
     +                                        conjg(acof(nb,lmp,nt)))
                                 ddnmt(llp,lh,nn) = ddnmt(llp,lh,nn) +
     +                                              we(nb)*real(cconst*
     +                                              bcof(nb,lm,nt)*
     +                                        conjg(bcof(nb,lmp,nt)))
                                 udnmt(llp,lh,nn) = udnmt(llp,lh,nn) +
     +                                              we(nb)*real(cconst*
     +                                              acof(nb,lm,nt)*
     +                                        conjg(bcof(nb,lmp,nt)))
                                 dunmt(llp,lh,nn) = dunmt(llp,lh,nn) +
     +                                              we(nb)*real(cconst*
     +                                              bcof(nb,lm,nt)*
     +                                        conjg(acof(nb,lmp,nt)))
   10                           ENDDO
                              ENDIF
   20                      ENDDO
   30                      natom = natom + neq(nn)
   40                   ENDDO
   50                ENDDO
                  ENDDO m_loop
   70          ENDDO
   80       ENDDO
   90    ENDDO
  100 ENDDO

      END SUBROUTINE rhonmt
      END MODULE m_rhonmt
