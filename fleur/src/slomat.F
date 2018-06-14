      MODULE m_slomat
c***********************************************************************
c updates the overlap matrix with the contributions from the local
c orbitals.
c                                                p.kurz sept. 1996
c***********************************************************************
      CONTAINS
      SUBROUTINE slomat(
     >                  llod,nlod,ntypd,natd,nvd,lmaxd,matsize,
     >                  ntyp,na,nv,invsat,fpi,con1,rmt,n_size,n_rank,
     >                  ddn,gk,rph,cph,f,g,kvec,
     >                  llo,nlo,dulon,uulon,uloulopn,
     >                  alo1,blo1,clo1,nlotot,l_ss,
     >                  ab_dim,iintsp,jintsp,chi11,chi22,chi21,
     X                  iilo,locol,nkvecprevat,bb)
c***********************************************************************
c iilo is the starting position for the LO-contribution of atom 'na'
c      in the Hamilton-matrix minus 1 (starting from (nv+1)*nv/2)
c nkvecprevat gives the number of G-vectors already processed on other
c             atoms (each l-LO gives (2l+1)*invsfact G-vectors)
c nkvecprevlo gives the number of G-vectors already processed on this
c             atoms due to other lo's (appears in hlomat only)
c locol stores the number of columns already processed; on parallel
c       computers this decides, whether the LO-contribution is
c       done on this node                                          gb00
!
! function legpol() at end of module
c***********************************************************************
#include"cpp_double.h"
c
      USE m_set, ONLY : cset
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nlod,ntypd,natd,nvd,lmaxd,matsize
      INTEGER, INTENT (IN) :: na,ntyp,n_size,n_rank,nlotot
      INTEGER, INTENT (IN) :: ab_dim,iintsp,jintsp
      LOGICAL, INTENT (IN) :: l_ss
      COMPLEX, INTENT (IN) :: chi11,chi22,chi21
      REAL,    INTENT (IN) :: con1,fpi
      INTEGER, INTENT (INOUT) :: iilo,nkvecprevat,locol
C     ..
C     .. Array Arguments ..
      INTEGER,INTENT (IN) :: invsat(natd),kvec(2* (2*llod+1),nlod)
      INTEGER,INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),nv(ab_dim)
      REAL,   INTENT (IN) :: alo1(nlod),blo1(nlod),clo1(nlod),rmt(ntypd)
      REAL,   INTENT (IN) :: ddn(0:lmaxd,ntypd),gk(nvd,3,ab_dim)
      REAL,   INTENT (IN) :: rph(nvd,ab_dim),cph(nvd,ab_dim)
      REAL,   INTENT (IN) :: f(nvd,0:lmaxd,ntypd,ab_dim)
      REAL,   INTENT (IN) :: g(nvd,0:lmaxd,ntypd,ab_dim)
      REAL,   INTENT (IN) :: dulon(nlod,ntypd),uulon(nlod,ntypd)
      REAL,   INTENT (IN) :: uloulopn(nlod,nlod,ntypd)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: bb(matsize)
#else
      COMPLEX, INTENT (INOUT) :: bb(matsize)
#endif
C     ..
C     .. Local Scalars ..
      REAL con,dotp,fact1,fact2,fact3,fl2p1
      INTEGER invsfct,k,kp,l,l2p1,lo,lop,lp,nkvec,nkvecat,nkvecp
      INTEGER iilo_s,locol_s,nkvecprevat_s,ic,ii,ij,n
      COMPLEX chihlp
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: bhelp(:)
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
c+noco
         IF (l_ss) THEN
           iilo_s = iilo
           locol_s = locol
           nkvecprevat_s = nkvecprevat ! save for other interstitial spin loops
           ic = 0                                        ! update b-matrix
           DO lo = 1,nlo(ntyp)
             ic = ic + invsfct* (2*llo(lo,ntyp)+1)
           ENDDO
           k = ic*(nv(jintsp)+nkvecprevat) + (ic+1)*ic/2
           ALLOCATE ( bhelp(k) )       ! initialize help-array 
           CALL cset(k,cmplx(0.,0.),bhelp,1)
           iilo = 0
         ENDIF
c-noco
         con = con1* ((rmt(ntyp))**2)/2.0
         nkvecat = 0
         DO lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            l2p1 = 2*l + 1
            fl2p1 = l2p1/fpi
#ifdef CPP_APW
            fact1 = (con**2)* fl2p1 * (
     +              alo1(lo)* (  alo1(lo) + 
     +                         2*clo1(lo)*   uulon(lo,ntyp) ) + 
     +              blo1(lo)* (  blo1(lo)*     ddn(l, ntyp) +
     +                         2*clo1(lo)*   dulon(lo,ntyp) ) +
     +              clo1(lo)*    clo1(lo)*uloulopn(lo,lo,ntyp) ) 
#else
            fact1 = (con**2)* fl2p1 * (
     +              alo1(lo)* (  alo1(lo) + 
     +                         2*clo1(lo) * uulon(lo,ntyp) ) +
     +              blo1(lo)* (  blo1(lo) *   ddn(l, ntyp) + 
     +                         2*clo1(lo) * dulon(lo,ntyp) ) +
     +              clo1(lo)*    clo1(lo) )
#endif
            DO nkvec = 1,invsfct* (2*l+1)
c+t3e
               locol = locol + 1
               IF (mod(locol-1,n_size).EQ.n_rank) THEN
c-t3e
               k = kvec(nkvec,lo)
c--->          calculate the overlap matrix elements with the regular
c--->          flapw basis-functions
               DO kp = 1,nv(jintsp)
                  iilo = iilo + 1
                  fact2 = con * fl2p1 * (
     +             f(kp,l,ntyp,jintsp)* ( alo1(lo) + 
     +                                    clo1(lo) * uulon(lo,ntyp) ) + 
     +             g(kp,l,ntyp,jintsp)* ( blo1(lo) *   ddn(l,ntyp) + 
     +                                    clo1(lo) * dulon(lo,ntyp) ) )
                  dotp = gk(k,1,iintsp) * gk(kp,1,jintsp) + 
     +                   gk(k,2,iintsp) * gk(kp,2,jintsp) +
     +                   gk(k,3,iintsp) * gk(kp,3,jintsp)
#ifdef CPP_INVERSION
                  bb(iilo) = bb(iilo) + invsfct*fact2 * legpol(l,dotp) *
     +                                  ( rph(k,iintsp)*rph(kp,jintsp) +
     +                                    cph(k,iintsp)*cph(kp,jintsp) )
#else
                  IF (.not.l_ss) THEN
                  bb(iilo) = bb(iilo) + invsfct*fact2 * legpol(l,dotp) *
     +                      cmplx( ( rph(k,iintsp)*rph(kp,jintsp) +
     +                               cph(k,iintsp)*cph(kp,jintsp) ) ,
     +                             ( cph(k,iintsp)*rph(kp,jintsp) -
     +                               rph(k,iintsp)*cph(kp,jintsp) ) )
                  ELSE 
                  bhelp(iilo) = invsfct*fact2*legpol(l,dotp) *
     +                      cmplx( ( rph(k,iintsp)*rph(kp,jintsp) +
     +                               cph(k,iintsp)*cph(kp,jintsp) ) ,
     +                             ( cph(k,iintsp)*rph(kp,jintsp) -
     +                               rph(k,iintsp)*cph(kp,jintsp) ) )
                  ENDIF
#endif
               END DO
c--->          calculate the overlap matrix elements with other local
c--->          orbitals at the same atom, if they have the same l
               iilo = iilo + nkvecprevat
               DO lop = 1, (lo-1)
                  lp = llo(lop,ntyp)
                  IF (l.EQ.lp) THEN
                     fact3 = con**2 * fl2p1 * (
     +                       alo1(lop)*(alo1(lo) + 
     +                                  clo1(lo)*   uulon(lo,ntyp)) +
     +                       blo1(lop)*(blo1(lo)*     ddn(l,ntyp) + 
     +                                  clo1(lo)*   dulon(lo,ntyp)) +
     +                       clo1(lop)*(alo1(lo)*   uulon(lop,ntyp) + 
     +                                  blo1(lo)*   dulon(lop,ntyp) + 
     +                                  clo1(lo)*uloulopn(lop,lo,ntyp)))
                     DO nkvecp = 1,invsfct* (2*lp+1)
                        iilo = iilo + 1
                        kp = kvec(nkvecp,lop) !
                        dotp = gk(k,1,iintsp) * gk(kp,1,jintsp) +
     +                         gk(k,2,iintsp) * gk(kp,2,jintsp) +
     +                         gk(k,3,iintsp) * gk(kp,3,jintsp)
#ifdef CPP_INVERSION
                        bb(iilo) =bb(iilo)+invsfct*fact3*legpol(l,dotp)* 
     +                                  ( rph(k,iintsp)*rph(kp,jintsp) +
     +                                    cph(k,iintsp)*cph(kp,jintsp) )
#else
                        IF (.not.l_ss) THEN
                        bb(iilo) =bb(iilo)+invsfct*fact3*legpol(l,dotp)*
     +                         cmplx( ( rph(k,iintsp)*rph(kp,jintsp) +
     +                                  cph(k,iintsp)*cph(kp,jintsp) ) ,
     +                                ( cph(k,iintsp)*rph(kp,jintsp) -
     +                                  rph(k,iintsp)*cph(kp,jintsp) ) )
                        ELSE 
                        bhelp(iilo) = invsfct*fact3*legpol(l,dotp)*
     +                         cmplx( ( rph(k,iintsp)*rph(kp,jintsp) +
     +                                  cph(k,iintsp)*cph(kp,jintsp) ) ,
     +                                ( cph(k,iintsp)*rph(kp,jintsp) -
     +                                  rph(k,iintsp)*cph(kp,jintsp) ) )
                        ENDIF
#endif
                     END DO
                  ELSE
                     iilo = iilo + invsfct* (2*lp+1)
                  END IF
               END DO
c--->          calculate the overlap matrix elements of one local
c--->          orbital with itself
               DO nkvecp = 1,nkvec
                  iilo = iilo + 1
                  kp = kvec(nkvecp,lo)
                  dotp = gk(k,1,iintsp) * gk(kp,1,jintsp) +
     +                   gk(k,2,iintsp) * gk(kp,2,jintsp) +
     +                   gk(k,3,iintsp) * gk(kp,3,jintsp)
#ifdef CPP_INVERSION
                  bb(iilo) = bb(iilo) + invsfct*fact1*legpol(l,dotp) *
     +                                ( rph(k,iintsp)*rph(kp,jintsp) +
     +                                  cph(k,iintsp)*cph(kp,jintsp) )
#else
                  IF (.not.l_ss) THEN
                  bb(iilo) = bb(iilo) + invsfct*fact1*legpol(l,dotp)*
     +                      cmplx( ( rph(k,iintsp)*rph(kp,jintsp) +
     +                               cph(k,iintsp)*cph(kp,jintsp) ) ,
     +                             ( cph(k,iintsp)*rph(kp,jintsp) -
     +                               rph(k,iintsp)*cph(kp,jintsp) ) )
                  ELSE 
                  bhelp(iilo) = invsfct*fact1*legpol(l,dotp)*
     +                      cmplx( ( rph(k,iintsp)*rph(kp,jintsp) +
     +                               cph(k,iintsp)*cph(kp,jintsp) ) ,
     +                             ( cph(k,iintsp)*rph(kp,jintsp) -
     +                               rph(k,iintsp)*cph(kp,jintsp) ) )
                  ENDIF
#endif
               END DO
               ENDIF ! mod(locol-1,n_size) = nrank 
c-t3e
            END DO
            nkvecat = nkvecat + invsfct* (2*l+1)
         END DO
         nkvecprevat = nkvecprevat + nkvecat
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
               CALL CPP_BLAS_caxpy(n-1,chihlp,bhelp(ii+1),1,bb(ij+1),1)
               ELSE
               CALL CPP_BLAS_caxpy(n,chihlp,bhelp(ii+1),1,bb(ij+1),1)
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
                 bb(ij) = bb(ij) +  chihlp *conjg( bhelp(ii) )
                 ij = ij + nv(1) + kp + nlotot
               ENDDO
             ENDDO
           ENDIF
c          n1=nv(jintsp)+1
c          n2=n1+nv(jintsp)+2
c          n3=n2+nv(jintsp)+3
c          do n = 1,nv(jintsp)+(nlotot/2)
c           write(4,'(8f15.8)') bhelp(n),bhelp(n1+n),bhelp(n2+n),bhelp(n3+n)
c          enddo
           DEALLOCATE ( bhelp )
           IF (.not.( iintsp.EQ.2 .AND. jintsp.EQ.2 )) THEN
             iilo = iilo_s
             locol = locol_s             ! restore for other loops
             nkvecprevat = nkvecprevat_s
           ENDIF
         ENDIF ! l_ss
c-noco
      END IF

      END SUBROUTINE slomat
!===========================================================================
      REAL FUNCTION legpol(l,arg)
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      REAL arg
      INTEGER l
C     ..
C     .. Local Scalars ..
      INTEGER lp
C     ..
C     .. Local Arrays ..
      REAL plegend(0:l)
C     ..
      plegend(0) = 1.0
      IF (l.GE.1) THEN
         plegend(1) = arg
         DO lp = 1,l - 1
            plegend(lp+1) = (lp+lp+1)*arg*plegend(lp)/ (lp+1) -
     +                      lp*plegend(lp-1)/ (lp+1)
         END DO
      END IF
      legpol = plegend(l)
      END FUNCTION legpol
!===========================================================================
      END MODULE m_slomat
