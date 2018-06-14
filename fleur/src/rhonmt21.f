      MODULE m_rhonmt21
c     *************************************************************
c     subroutine sets up the coefficients of the spin (up,down) 
c     part of the non-spherical muffin-tin density. 
c                                                 pk`00 ff`01 gb`02
c     *************************************************************
      CONTAINS
      SUBROUTINE rhonmt21(
     >                    lmaxd,llpd,ntypd,memd,nlhd,ntypsd,natd,
     >                    nobd,jspd,we,ne,ntype,nsymt,neq,lmax,
     >                    clnu,nmem,mlh,nlh,llh,ntypsy,
     >                    acof,bcof,
     X                    uunmt21,ddnmt21,udnmt21,dunmt21)
      USE m_gaunt,ONLY:gaunt1
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT(IN) :: lmaxd,llpd,ntypd,memd,nlhd,ntypsd,natd,nobd
      INTEGER,INTENT(IN) :: ne,ntype,nsymt,jspd
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT(IN) :: neq(ntypd),ntypsy(natd),lmax(ntypd)
      INTEGER, INTENT(IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT(IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      COMPLEX, INTENT(IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT(IN) :: acof(nobd,0:lmaxd* (lmaxd+2),natd,jspd)
      COMPLEX, INTENT(IN) :: bcof(nobd,0:lmaxd* (lmaxd+2),natd,jspd)
      REAL,    INTENT(IN) :: we(nobd)
      COMPLEX, INTENT (INOUT) :: ddnmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX, INTENT (INOUT) :: dunmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX, INTENT (INOUT) :: udnmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX, INTENT (INOUT) :: uunmt21((lmaxd+1)**2,nlhd,ntypd)
C     ..
C     .. Local Scalars ..
      COMPLEX coef, cconst, cil, coef1
      COMPLEX, PARAMETER :: mi = cmplx(0.0,-1.0)
      INTEGER jmem,l,lh,llp,lm,lmp,lp,lv,
     +        m,mp,mv,na,natom,nb,nn,ns,nt
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,mod,real,cmplx
C     ..
c
      DO ns=1,nsymt
        natom= 0
        DO nn=1,ntype
          nt= natom
          DO na= 1,neq(nn)
            nt= nt+1
            IF (ntypsy(nt)==ns) THEN

              DO lh = 1,nlh(ns)
                lv = llh(lh,ns)
                DO lp = 0,lmax(nn)
                  DO l = 0,lmax(nn)

                    IF ( MOD(lv+l+lp,2) == 0 ) THEN
                      cil = mi**(l-lp)
                      llp= lp*(lmax(nn)+1)+l+1

                      DO jmem = 1,nmem(lh,ns)
                        mv = mlh(jmem,lh,ns)
                        coef1 = cil * clnu(jmem,lh,ns) 
                        DO mp = -lp,lp
                          lmp = lp*(lp+1) + mp
                          DO m = -l,l
                            lm= l*(l+1) + m
                            coef=  conjg( coef1 *               
     *                             gaunt1(l,lv,lp,m,mv,mp,lmaxd) )

                            IF (abs(coef) >= 0 ) THEN
                              DO nb = 1,ne
                                cconst= we(nb) * coef
                                uunmt21(llp,lh,nn) = uunmt21(llp,lh,nn)+ 
     +                                        cconst * acof(nb,lm,nt,1)*
     *                                          conjg(acof(nb,lmp,nt,2))
                                udnmt21(llp,lh,nn) = udnmt21(llp,lh,nn)+
     +                                        cconst * bcof(nb,lm,nt,1)*
     *                                          conjg(acof(nb,lmp,nt,2))
                                dunmt21(llp,lh,nn) = dunmt21(llp,lh,nn)+
     +                                        cconst * acof(nb,lm,nt,1)*
     *                                          conjg(bcof(nb,lmp,nt,2))
                                ddnmt21(llp,lh,nn) = ddnmt21(llp,lh,nn)+
     +                                        cconst * bcof(nb,lm,nt,1)*
     *                                          conjg(bcof(nb,lmp,nt,2))
                              ENDDO ! nb
                            ENDIF ! (coef >= 0)

                          ENDDO ! mp
                        ENDDO ! m
                      ENDDO ! jmem

                    ENDIF ! ( MOD(lv+l+lp),2) == 0 )

                  ENDDO ! lp
                ENDDO ! l
              ENDDO ! lh

            ENDIF ! (ntypsy(nt)==ns)
          ENDDO ! na
          natom= natom + neq(nn)
        ENDDO ! nn

      ENDDO ! ns

      RETURN
      
      END SUBROUTINE rhonmt21
      END MODULE m_rhonmt21
