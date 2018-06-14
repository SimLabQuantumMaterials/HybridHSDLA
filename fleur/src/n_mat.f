      MODULE m_nmat
c     ************************************************************
c     This subroutine calculates the density matrix n^{s}_{m,m'}
c     for a given atom 'n' and l-quantum number 'l'. The l's for
c     all atoms are stored in lda_u(), if lda_u()<0, no +U is used.
c     For details see Eq.(12) of Shick et al. PRB 60, 10765 (1999)
c     Part of the LDA+U package                   G.B., Oct. 2000
c     ************************************************************
      CONTAINS
      SUBROUTINE n_mat(
     >                 lmaxd,ntypd,neigd,nobd,natd,nop,llod,nlod,
     >                 ngopr,neq,lmax,ne,ddn,we,ntype,lda_u,n_u,
     >                 uulon,dulon,llo,nlo,uloulopn,invsat,invtab,
     >                 lmd,d_wgn,acof,bcof,ccof,invarind,invarop,
     <                 n_mmp)
C
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,neigd,nobd,natd
      INTEGER, INTENT (IN) :: llod,nlod,nop,lmd
      INTEGER, INTENT (IN) :: ntype,ne,n_u
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ngopr(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: neq(ntypd),lda_u(ntype)
      INTEGER, INTENT (IN) :: invsat(natd),invtab(nop)
      INTEGER, INTENT (IN) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: uulon(nlod,ntypd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd),we(neigd)
      COMPLEX, INTENT (IN) :: d_wgn(-3:3,-3:3,3,nop)
      COMPLEX, INTENT (IN) :: acof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (INOUT) :: n_mmp(-3:3,-3:3,n_u)
C     ..
C     .. Local Scalars ..
      COMPLEX c_0
      INTEGER i,j,k,l,m,mp,n,it,is,isi,natom,n_ldau,lp
      INTEGER ilo,ilop,ll1,nn,lmp,lm
      REAL fac
C     ..
C     .. Local Arrays ..
      COMPLEX n_tmp(-3:3,-3:3),nr_tmp(-3:3,-3:3),d_tmp(-3:3,-3:3)
      COMPLEX n1_tmp(-3:3,-3:3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg
C     ..
c
c calculate n_mat:
c
      n_ldau = 0
      natom = 0
      DO n = 1,ntype
        IF (lda_u(n).GE.0) THEN
          n_ldau = n_ldau + 1
          n_tmp(:,:) =cmplx(0.0,0.0)
          l = lda_u(n)
          ll1 = (l+1)*l 
          DO nn = 1, neq(n)
            natom = natom + 1
!
!  prepare n_mat in local frame (in noco-calculations this depends 
!                                also on alpha(n) and beta(n) )
!
            DO m = -l,l
              lm = ll1+m
              DO mp = -l,l
                lmp = ll1+mp
                c_0 = cmplx(0.0,0.0)
                DO i = 1,ne
                  c_0 = c_0 +  we(i) * ( ddn(l,n) *
     +              conjg(bcof(i,lmp,natom))*bcof(i,lm,natom) +
     +              conjg(acof(i,lmp,natom))*acof(i,lm,natom) )
                ENDDO
                n_tmp(m,mp) = c_0 
              ENDDO
            ENDDO
!
!  add local orbrbital contribution (if there is one) (untested so far)
!
            DO ilo = 1, nlo(n)
              IF (llo(ilo,n).EQ.l) THEN

                 DO m = -l,l
                   lm = ll1+m
                   DO mp = -l,l
                     lmp = ll1+mp
                     c_0 = cmplx(0.0,0.0)
                     DO i = 1,ne
                       c_0 = c_0 +  we(i) * (  uulon(ilo,n) * (
     +                    conjg(acof(i,lmp,natom))*ccof(m,i,ilo,natom) +
     +                    conjg(ccof(mp,i,ilo,natom))*acof(i,lm,natom) )
     +                                       + dulon(ilo,n) * (
     +                    conjg(bcof(i,lmp,natom))*ccof(m,i,ilo,natom) +
     +                    conjg(ccof(mp,i,ilo,natom))*bcof(i,lm,natom)))
                     ENDDO
                     DO ilop = 1, nlo(n)
                       IF (llo(ilop,n).EQ.l) THEN
                         DO i = 1,ne
                           c_0 = c_0 +  we(i) * uloulopn(ilo,ilop,n) *
     +                                  conjg(ccof(mp,i,ilop,natom)) *
     +                                        ccof(m ,i,ilo ,natom)
                         ENDDO 
                       ENDIF
                     ENDDO
                     n_tmp(m,mp) = n_tmp(m,mp) + c_0
                   ENDDO
                 ENDDO

              ENDIF
            ENDDO
!
!  n_mmp should be rotated by D_mm' ; compare force_a21
!
            DO it = 1, invarind(natom)

              fac = 1.0  /  ( invarind(natom) * neq(n) )
               is = invarop(natom,it)
              isi = invtab(is)
              d_tmp(:,:) = cmplx(0.0,0.0)
              DO m = -l,l
                DO mp = -l,l
                  d_tmp(m,mp) = d_wgn(m,mp,l,isi)
                ENDDO
              ENDDO
              nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
              n1_tmp =  matmul( nr_tmp, d_tmp )
              DO m = -l,l
                DO mp = -l,l
                  n_mmp(m,mp,n_ldau) = n_mmp(m,mp,n_ldau) +
     +                                conjg(n1_tmp(m,mp)) * fac
                ENDDO
              ENDDO

            ENDDO

          ENDDO ! sum  over equivalent atoms
        ELSE
          natom = natom + neq(n)
        ENDIF
      ENDDO     ! loop over atom types
    
!     do m=-l,l
!      write(*,'(14f12.6)') (n_mmp(m,mp),mp=-l,l)
!     enddo
c
      RETURN
      END SUBROUTINE n_mat
      END MODULE m_nmat
