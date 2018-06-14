      MODULE m_hlomat
c***********************************************************************
c updates the hamiltonian  matrix with the contributions from the local
c orbitals.
c p.kurz sept. 1996
c***********************************************************************
      CONTAINS
      SUBROUTINE hlomat(
     >                  llod,nlod,ntypd,natd,nvd,
     >                  matsize,loplod,lmd,lmplmd,n_size,n_rank,
     >                  ntyp,na,nv,lmax,invsat,ar,br,ai,bi,
     >                  tuu,tud,tdu,tdd,ind,el,ddn,
     >                  llo,nlo,alo,blo,clo,tuulo,tdulo,tuloulo,
     >                  lmaxd,us,dus,uds,duds,ulos,dulos,rmt,nlotot,
     >                  l_ss,ab_dim,iintsp,jintsp,chi11,chi22,chi21,
     X                  iilo,locol,nkvecprevat,aa)
#include"cpp_double.h"
c
      USE m_set, ONLY : cset
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nlod,ntypd,natd,nvd,matsize
      INTEGER, INTENT (IN) :: loplod,lmd,lmplmd,n_size,n_rank
      INTEGER, INTENT (IN) :: lmax,na,ntyp,lmaxd,nlotot
      INTEGER, INTENT (IN) :: ab_dim,iintsp,jintsp
      LOGICAL, INTENT (IN) :: l_ss
      COMPLEX, INTENT (IN) :: chi11,chi22,chi21
      INTEGER, INTENT (INOUT) :: iilo,nkvecprevat,locol
      REAL,    INTENT (IN) :: rmt
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd)
      INTEGER, INTENT (INOUT) :: ind(0:lmd,0:lmd)
      INTEGER, INTENT (IN) :: invsat(natd),nv(ab_dim)
      COMPLEX, INTENT (IN) :: tdulo(0:lmd,-llod:llod,nlod)
      COMPLEX, INTENT (IN) :: tuloulo(-llod:llod,-llod:llod,loplod)
      COMPLEX, INTENT (IN) :: tuulo(0:lmd,-llod:llod,nlod)
      COMPLEX, INTENT (IN) :: alo(-llod:llod,2* (2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT (IN) :: blo(-llod:llod,2* (2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT (IN) :: clo(-llod:llod,2* (2*llod+1),nlod,ab_dim)
      REAL,    INTENT (IN) :: ar(nvd,0:lmd),br(nvd,0:lmd)
      REAL,    INTENT (IN) :: ai(nvd,0:lmd),bi(nvd,0:lmd)
      COMPLEX, INTENT (INOUT) :: tdd(0:lmplmd),tdu(0:lmplmd)
      COMPLEX, INTENT (INOUT) :: tud(0:lmplmd),tuu(0:lmplmd)
      REAL,    INTENT (IN) :: el(0:lmax),ddn(0:lmax)
      REAL,    INTENT (IN) :: us(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: dus(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: ulos(nlod,ntypd),dulos(nlod,ntypd)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: aa(matsize)
#else
      COMPLEX, INTENT (INOUT) :: aa(matsize)
#endif
C     ..
C     .. Local Scalars ..
      COMPLEX axx,bxx,cxx,dtd,dtu,dtulo,ulotd,ulotu,ulotulo,utd,utu,
     +        utulo,chihlp
      INTEGER im,in,invsfct,kp,l,lm,lmp,lo,lolo,lolop,lop,lp,m,matel0,
     +        locol0,mp,nkvec,nkvecp,nkvecprevlo,lmplm 
      INTEGER iilo_s,locol_s,nkvecprevat_s,ic,ii,ij,n,k
C     ..
C     .. Local Arrays ..
      COMPLEX ax(nvd),bx(nvd),cx(nvd)
      COMPLEX, ALLOCATABLE :: ahelp(:)
      INTEGER indt(0:lmd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aimag,cmplx,conjg,real
C     ..
C     .. External Subroutines ..
      EXTERNAL CPP_BLAS_caxpy
C     ..

      IF ((invsat(na).EQ.0) .OR. (invsat(na).EQ.1)) THEN
c--->    if this atom is the first of two atoms related by inversion,
c--->    the contributions to the overlap matrix of both atoms are added
c--->    at once. where it is made use of the fact, that the sum of
c--->    these contributions is twice the real part of the contribution
c--->    of each atom. note, that in this case there are twice as many
c--->    (2*(2*l+1)) k-vectors (compare abccoflo and comments there).
         IF (invsat(na).EQ.0) invsfct = 1
         IF (invsat(na).EQ.1) invsfct = 2
         nkvecprevlo = 0
c-noco
         IF (l_ss) THEN
           iilo_s = iilo
           locol_s = locol
           nkvecprevat_s = nkvecprevat ! save for other interstitial spin loops
           ic = 0                                        ! update b-matrix
           DO lo = 1,nlo(ntyp)
             ic = ic + invsfct* (2*llo(lo,ntyp)+1)
           ENDDO
           k = ic*(nv(jintsp)+nkvecprevat) + (ic+1)*ic/2
           ALLOCATE ( ahelp(k) )       ! initialize help-array
           CALL cset(k,cmplx(0.,0.),ahelp,1)
           iilo = 0
         ENDIF
c-noco
c
c temporarily update the diagonal elements
c
         DO l = 0,lmax
            DO  m = -l,l
              lm = l* (l+1) + m
              lmplm = (lm* (lm+3))/2
              tuu(lmplm) = tuu(lmplm) + el(l)
              tdd(lmplm) = tdd(lmplm) + el(l)*ddn(l)
              tud(lmplm) = tud(lmplm) + 0.5
              tdu(lmplm) = tdu(lmplm) + 0.5
              indt(lm) = ind(lm,lm)
              ind(lm,lm) = lmplm
            ENDDO
         ENDDO
c
         DO lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            matel0 = iilo
            locol0 = locol
c--->       calculate the hamiltonian matrix elements with the regular
c--->       flapw basis-functions
            DO m = -l,l
               lm = l* (l+1) + m
               iilo = matel0
               locol = locol0
               DO kp = 1,nv(jintsp)
                  ax(kp) = cmplx(0.0,0.0)
                  bx(kp) = cmplx(0.0,0.0)
                  cx(kp) = cmplx(0.0,0.0)
               END DO
               DO lp = 0,lmax
                  DO mp = -lp,lp
                     lmp = lp* (lp+1) + mp
                     in = ind(lmp,lm)
                     IF (in.NE.-9999) THEN
                        IF (in.GE.0) THEN
                           utu = tuu(in)
                           dtu = tdu(in)
                           utd = tud(in)
                           dtd = tdd(in)
                        ELSE
                           im = -in
                           utu = conjg(tuu(im))
                           dtu = conjg(tud(im))
                           utd = conjg(tdu(im))
                           dtd = conjg(tdd(im))
                        END IF
                        utulo = tuulo(lmp,m,lo)
                        dtulo = tdulo(lmp,m,lo)
c--->                   note, that utu,dtu... are the t-matrices and
c--->                   not their complex conjugates as in hssphn
c--->                   and that a,b,alo... are the complex
c--->                   conjugates of the a,b...-coefficients
                        DO kp = 1,nv(jintsp)
                           ax(kp) = ax(kp) + 
     +                              cmplx(ar(kp,lmp),ai(kp,lmp))*utu +
     +                              cmplx(br(kp,lmp),bi(kp,lmp))*dtu
                           bx(kp) = bx(kp) + 
     +                              cmplx(ar(kp,lmp),ai(kp,lmp))*utd +
     +                              cmplx(br(kp,lmp),bi(kp,lmp))*dtd
                           cx(kp) = cx(kp) + 
     +                              cmplx(ar(kp,lmp),ai(kp,lmp))*utulo +
     +                              cmplx(br(kp,lmp),bi(kp,lmp))*dtulo
                        END DO
                     END IF
                  END DO
               END DO
c+t3e
               DO nkvec = 1,invsfct* (2*l+1)
                  locol = locol + 1
                  IF (mod(locol-1,n_size).EQ.n_rank) THEN
c-t3e
                  DO kp = 1,nv(jintsp)
                     iilo = iilo + 1
#ifdef CPP_INVERSION
                     aa(iilo) = aa(iilo) + invsfct * (
     +                      real(alo(m,nkvec,lo,iintsp))* real(ax(kp)) +
     +                     aimag(alo(m,nkvec,lo,iintsp))*aimag(ax(kp)) +
     +                      real(blo(m,nkvec,lo,iintsp))* real(bx(kp)) +
     +                     aimag(blo(m,nkvec,lo,iintsp))*aimag(bx(kp)) +
     +                      real(clo(m,nkvec,lo,iintsp))* real(cx(kp)) +
     +                     aimag(clo(m,nkvec,lo,iintsp))*aimag(cx(kp)) )
c---> APWlo
                     aa(iilo) = aa(iilo) + 0.25 * rmt**2 * invsfct * (
     +               (ar(kp,lm)*us(l,ntyp) + br(kp,lm)*uds(l,ntyp))*
     +                   real( alo(m,nkvec,lo,iintsp)*  dus(l,ntyp)  +
     +                         blo(m,nkvec,lo,iintsp)* duds(l,ntyp)  +
     +                         clo(m,nkvec,lo,iintsp)*dulos(lo,ntyp) ) +
     +               (ai(kp,lm)*us(l,ntyp) + bi(kp,lm)*uds(l,ntyp))*
     +                  aimag( alo(m,nkvec,lo,iintsp)*  dus(l,ntyp) +
     +                         blo(m,nkvec,lo,iintsp)* duds(l,ntyp) +
     +                         clo(m,nkvec,lo,iintsp)*dulos(lo,ntyp) ) )
#else
                     IF (.not.l_ss) THEN
                     aa(iilo) = aa(iilo) + invsfct * (
     +                        alo(m,nkvec,lo,iintsp) * conjg( ax(kp) ) +
     +                        blo(m,nkvec,lo,iintsp) * conjg( bx(kp) ) +
     +                        clo(m,nkvec,lo,iintsp) * conjg( cx(kp) ) )
c---> APWlo
                     aa(iilo)=aa(iilo) + 0.25 * rmt**2 * invsfct*(
     $                      (cmplx(ar(kp,lm),-ai(kp,lm))* us(l,ntyp)+
     +                       cmplx(br(kp,lm),-bi(kp,lm))*uds(l,ntyp))*
     $                       (alo(m,nkvec,lo,iintsp)*  dus(l,ntyp)
     +                       +blo(m,nkvec,lo,iintsp)* duds(l,ntyp)
     +                       +clo(m,nkvec,lo,iintsp)*dulos(lo,ntyp) ) )
                     ELSE
                     ahelp(iilo) = ahelp(iilo) + invsfct * (
     +                        alo(m,nkvec,lo,iintsp) * conjg( ax(kp) ) +
     +                        blo(m,nkvec,lo,iintsp) * conjg( bx(kp) ) +
     +                        clo(m,nkvec,lo,iintsp) * conjg( cx(kp) ) )
c---> APWlo
                     ahelp(iilo)=ahelp(iilo) + 0.25 * rmt**2 * invsfct*(
     $                      (cmplx(ar(kp,lm),-ai(kp,lm))* us(l,ntyp)+
     +                       cmplx(br(kp,lm),-bi(kp,lm))*uds(l,ntyp))*
     $                       (alo(m,nkvec,lo,iintsp)*  dus(l,ntyp)
     +                       +blo(m,nkvec,lo,iintsp)* duds(l,ntyp)
     +                       +clo(m,nkvec,lo,iintsp)*dulos(lo,ntyp) ) )
                     ENDIF
#endif
                  END DO
c--->             jump to the last matrixelement of the current row
                  iilo = iilo + nkvecprevat + nkvecprevlo + nkvec
                  ENDIF
               END DO
            END DO
c--->       calculate the hamiltonian matrix elements with other
c--->       local orbitals at the same atom and with itself
            iilo = matel0
c+t3e
            locol = locol0
            DO nkvec = 1,invsfct* (2*l+1)
               locol = locol + 1
               IF (mod(locol-1,n_size).EQ.n_rank) THEN
c-t3e
c--->          skip the matrixelements with regular flapw-fcn. and
c--->          with local orbitals at other atoms
               iilo = iilo + nv(jintsp) + nkvecprevat
c--->          calculate the hamiltonian matrix elements with other
c--->          local orbitals at the same atom, if they have the same l
               DO lop = 1, (lo-1)
                  lp = llo(lop,ntyp)
                  DO nkvecp = 1,invsfct* (2*lp+1)
                     iilo = iilo + 1
                     DO m = -l,l
                        lm = l* (l+1) + m
                        DO mp = -lp,lp
                           lmp = lp* (lp+1) + mp
                           in = ind(lmp,lm)
                           IF (in.NE.-9999) THEN
                              IF (in.GE.0) THEN
                                utu = tuu(in)
                                dtu = tdu(in)
                                utd = tud(in)
                                dtd = tdd(in)
                              ELSE
                                im = -in
                                utu = conjg(tuu(im))
                                dtu = conjg(tud(im))
                                utd = conjg(tdu(im))
                                dtd = conjg(tdd(im))
                              END IF
                              utulo = tuulo(lmp,m,lo)
                              dtulo = tdulo(lmp,m,lo)
                              ulotu = conjg(tuulo(lm,mp,lop))
                              ulotd = conjg(tdulo(lm,mp,lop))
c--->                         note that lo > lop
                              lolop = (lo-1)*lo/2 + lop
                              ulotulo = conjg(tuloulo(m,mp,lolop))
                              axx=conjg(alo(m,nkvec,lo,iintsp))*utu +
     +                            conjg(blo(m,nkvec,lo,iintsp))*utd +
     +                            conjg(clo(m,nkvec,lo,iintsp))*utulo
                              bxx=conjg(alo(m,nkvec,lo,iintsp))*dtu +
     +                            conjg(blo(m,nkvec,lo,iintsp))*dtd +
     +                            conjg(clo(m,nkvec,lo,iintsp))*dtulo
                              cxx = 
     +                          conjg(alo(m,nkvec,lo,iintsp))*ulotu +
     +                          conjg(blo(m,nkvec,lo,iintsp))*ulotd +
     +                          conjg(clo(m,nkvec,lo,iintsp))*ulotulo
#ifdef CPP_INVERSION
                                 aa(iilo) = aa(iilo) + invsfct * (
     +                      real(alo(mp,nkvecp,lop,jintsp))* real(axx) -
     +                     aimag(alo(mp,nkvecp,lop,jintsp))*aimag(axx) +
     +                      real(blo(mp,nkvecp,lop,jintsp))* real(bxx) -
     +                     aimag(blo(mp,nkvecp,lop,jintsp))*aimag(bxx) +
     +                      real(clo(mp,nkvecp,lop,jintsp))* real(cxx) -
     +                     aimag(clo(mp,nkvecp,lop,jintsp))*aimag(cxx) )
#else
                              IF (.not.l_ss) THEN
                                 aa(iilo) = aa(iilo) + invsfct * conjg(
     +                                 alo(mp,nkvecp,lop,jintsp) * axx +
     +                                 blo(mp,nkvecp,lop,jintsp) * bxx +
     +                                 clo(mp,nkvecp,lop,jintsp) * cxx )
                              ELSE
                                 ahelp(iilo)=ahelp(iilo)+invsfct*conjg(
     +                                 alo(mp,nkvecp,lop,jintsp) * axx +
     +                                 blo(mp,nkvecp,lop,jintsp) * bxx +
     +                                 clo(mp,nkvecp,lop,jintsp) * cxx )
                              ENDIF
#endif
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
c--->          calculate the hamiltonian matrix elements of one local
c--->          orbital with itself
               DO nkvecp = 1,nkvec
                  iilo = iilo + 1
                  DO m = -l,l
                     lm = l* (l+1) + m
                     DO mp = -l,l
                        lmp = l* (l+1) + mp
                        in = ind(lmp,lm)
                        IF (in.NE.-9999) THEN
                           IF (in.GE.0) THEN
                              utu = tuu(in)
                              dtu = tdu(in)
                              utd = tud(in)
                              dtd = tdd(in)
                            ELSE
                              im = -in
                              utu = conjg(tuu(im))
                              dtu = conjg(tud(im))
                              utd = conjg(tdu(im))
                              dtd = conjg(tdd(im))
                           END IF
                           utulo = tuulo(lmp,m,lo)
                           dtulo = tdulo(lmp,m,lo)
                           ulotu = conjg(tuulo(lm,mp,lo))
                           ulotd = conjg(tdulo(lm,mp,lo))
                           lolo = (lo-1)*lo/2 + lo
                           ulotulo = conjg(tuloulo(m,mp,lolo))
                           axx = conjg(alo(m,nkvec,lo,iintsp))*utu +
     +                           conjg(blo(m,nkvec,lo,iintsp))*utd +
     +                           conjg(clo(m,nkvec,lo,iintsp))*utulo
                           bxx = conjg(alo(m,nkvec,lo,iintsp))*dtu +
     +                           conjg(blo(m,nkvec,lo,iintsp))*dtd +
     +                           conjg(clo(m,nkvec,lo,iintsp))*dtulo
                           cxx = conjg(alo(m,nkvec,lo,iintsp))*ulotu +
     +                           conjg(blo(m,nkvec,lo,iintsp))*ulotd +
     +                           conjg(clo(m,nkvec,lo,iintsp))*ulotulo
#ifdef CPP_INVERSION
                           aa(iilo) = aa(iilo) + invsfct* (
     +                       real(alo(mp,nkvecp,lo,jintsp))* real(axx) -
     +                      aimag(alo(mp,nkvecp,lo,jintsp))*aimag(axx) +
     +                       real(blo(mp,nkvecp,lo,jintsp))* real(bxx) -
     +                      aimag(blo(mp,nkvecp,lo,jintsp))*aimag(bxx) +
     +                       real(clo(mp,nkvecp,lo,jintsp))* real(cxx) -
     +                      aimag(clo(mp,nkvecp,lo,jintsp))*aimag(cxx) )
#else
                           IF (.not.l_ss) THEN
                           aa(iilo) = aa(iilo) + invsfct* conjg(
     +                                    alo(mp,nkvecp,lo,jintsp)*axx +
     +                                    blo(mp,nkvecp,lo,jintsp)*bxx +
     +                                    clo(mp,nkvecp,lo,jintsp)*cxx )
                           ELSE
                           ahelp(iilo) = ahelp(iilo) + invsfct* conjg(
     +                                    alo(mp,nkvecp,lo,jintsp)*axx +
     +                                    blo(mp,nkvecp,lo,jintsp)*bxx +
     +                                    clo(mp,nkvecp,lo,jintsp)*cxx )
                           ENDIF
#endif
                        END IF
                     END DO
                  END DO
               END DO
               ENDIF
c-t3e
            END DO
            nkvecprevlo = nkvecprevlo + invsfct* (2*l+1)
         END DO ! end of lo = 1,nlo loop
c
c remove the temporary update of the diagonal elements
c
         DO l = 0,lmax
            DO  m = -l,l
              lm = l* (l+1) + m
              lmplm = (lm* (lm+3))/2
              tuu(lmplm) = tuu(lmplm) - el(l)
              tdd(lmplm) = tdd(lmplm) - el(l)*ddn(l)
              tud(lmplm) = tud(lmplm) - 0.5
              tdu(lmplm) = tdu(lmplm) - 0.5
              ind(lm,lm) = indt(lm)
            ENDDO
         ENDDO
c
         nkvecprevat = nkvecprevat + nkvecprevlo
c+noco
         IF (l_ss) THEN
           IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN     !---> spin-up spin-up part
              chihlp = chi11
              n = nv(1) + nkvecprevat_s
              ij = (n+1)*n/2
           ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN !---> spin-down spin-down part
              chihlp = chi22
              n = nv(1) + nv(iintsp) + nlotot + nkvecprevat_s
              ij = (n+1)*n/2 + nv(1) + nlotot
           ELSE                                          !---> spin-down spin-up part
              chihlp = chi21
              n = nv(1) + nv(iintsp) + nlotot + nkvecprevat_s
              ij = (n+1)*n/2
           ENDIF

           ic = 0                                        ! update b-matrix
           DO lo = 1,nlo(ntyp)
             ic = ic + invsfct* (2*llo(lo,ntyp)+1)
           ENDDO
           IF (.not.( iintsp.EQ.1 .AND. jintsp.EQ.2 )) THEN
             ii = 0
             DO k = 1, ic
               n = k + nv(jintsp) + nkvecprevat_s
#ifndef CPP_INVERSION
               IF (iintsp.EQ.2 .AND. jintsp.EQ.1 ) THEN
               CALL CPP_BLAS_caxpy(n-1,chihlp,ahelp(ii+1),1,aa(ij+1),1)
               ELSE
               CALL CPP_BLAS_caxpy(n,chihlp,ahelp(ii+1),1,aa(ij+1),1)
               ENDIF
#endif
               ii = ii + n
               ij = ij + n + (nv(3-jintsp)+nlotot)*(iintsp-1)
             ENDDO
           ELSE                                         ! special treatment for up-down:
             n = nv(1) + nlotot
             ii = 0
             DO k = 1, ic
               ij = (n+1)*n/2 + nv(1) + k + nkvecprevat_s
               DO kp = 1, k + nv(jintsp) + nkvecprevat_s
                 ii = ii + 1
                 aa(ij) = aa(ij) + chihlp * conjg( ahelp(ii) )
                 ij = ij + nv(1) + kp + nlotot
               ENDDO
             ENDDO
           ENDIF
           DEALLOCATE ( ahelp )
           IF (.not.( iintsp.EQ.2 .AND. jintsp.EQ.2 )) THEN
             iilo = iilo_s
             locol = locol_s             ! restore for other loops
             nkvecprevat = nkvecprevat_s
           ENDIF
         ELSE
           k = nv(1) * ( nv(1) + 1 ) / 2
         ENDIF
c-noco
      END IF
      END SUBROUTINE hlomat
      END MODULE m_hlomat
