      MODULE m_uham
c-------------------------------------------------------------------+
c     For details see Eq.(21) of Shick et al. PRB 60, 10765 (1999)  |
c     Sets up the LDA + U Hamilton matrix with eigenvalues in the   |
c     diagonal:                                                     |
c  s     ---       *s,j'  s,j    . .   *s,j'  s,j    s            s |
c H    = >      [ A      A    + <u|u> B      B    ] V     + d    e  |
c  j',j  -- m,m'   l,m'   l,m          l,m'   l,m    m,m'    j,j' j |
c                                                                   |
C                                                  G.B. Oct. 2000   |
c                                                                   |
c Version for 1st variation to be called from hssphn.               |
c now the Hamiltonmatrix is close packed aa()  and ar,ai and br,bi  |
c are provided instead of ahelp and bhelp.              Feb. 2001   |
c-------------------------------------------------------------------+
      CONTAINS
      SUBROUTINE u_ham(
     >                 jspd,nvd,lmaxd,ntypd,lmd,matsize,llod,
     >                 nlod,jspins,nv,lda_u,jsp,n,invsfct,
     >                 ar,ai,br,bi,ddn,vs_mmp,n_u,lmaxb,
     >                 llo,nlo,alo,blo,clo,uulon,dulon,uloulopn,
     >                 n_size,n_rank,nlotot,
     >                 l_ss,ab_dim,iintsp,jintsp,chi11,chi22,chi21,
     X                 nkvecprevat,iilo,locol,aa)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,nvd,lmaxd,ntypd,lmaxb,lmd,llod
      INTEGER, INTENT (IN) :: nlod,jspins,n_u,jsp,n,matsize
      INTEGER, INTENT (IN) :: n_size,n_rank,nlotot
      INTEGER, INTENT (IN) :: ab_dim,iintsp,jintsp
      INTEGER, INTENT (INOUT) :: nkvecprevat,iilo,locol
      REAL,    INTENT (IN) :: invsfct
      LOGICAL, INTENT (IN) :: l_ss
      COMPLEX, INTENT (IN) :: chi11,chi22,chi21
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lda_u(ntypd),nv(jspd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN) :: ar(nvd,0:lmd,ab_dim),ai(nvd,0:lmd,ab_dim)
      REAL,    INTENT (IN) :: br(nvd,0:lmd,ab_dim),bi(nvd,0:lmd,ab_dim)
      COMPLEX, INTENT (IN) :: alo(-llod:llod,2* (2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT (IN) :: blo(-llod:llod,2* (2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT (IN) :: clo(-llod:llod,2* (2*llod+1),nlod,ab_dim)
      COMPLEX,INTENT(IN):: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: aa(matsize)
#else
      COMPLEX, INTENT (INOUT) :: aa(matsize)
#endif
C     ..
C     .. Local Scalars ..
      INTEGER m,mp,ispin,itype,l,ig,igp,ll,llm,nvv,ii,n_l,lp,lop
      INTEGER matel0,locol0,lo,kp,nkvec,nkvecp,nkvecprevlo,nkend
      INTEGER igp_m,i_invsf,nc
      COMPLEX chihlp
      REAL    dd,uulo,dulo,ulolo
C     ..
C     .. Local Arrays
      REAL, ALLOCATABLE :: ar_help(:,:),ai_help(:,:),va_help(:)
      REAL, ALLOCATABLE :: br_help(:,:),bi_help(:,:),vb_help(:)
      REAL, ALLOCATABLE :: cr_help(:,:),ci_help(:,:)
      COMPLEX, ALLOCATABLE :: ax(:,:),bx(:,:),cx(:,:)
      COMPLEX, ALLOCATABLE :: wa_help(:),wb_help(:),wc_help(:)
      COMPLEX, ALLOCATABLE :: wd_help(:),we_help(:)


      INTRINSIC cmplx,conjg

      i_invsf = nint(invsfct)
      n_l = 0 
      DO itype = 1,n
        l = lda_u(itype)
        IF (l.GE.0) n_l = n_l + 1
      ENDDO
      nvv = nv(jsp)
      IF (l_ss)  nvv = nv(jintsp)
      
      ALLOCATE ( ar_help(nvv,-l:l),ai_help(nvv,-l:l),va_help(-l:l),
     +           br_help(nvv,-l:l),bi_help(nvv,-l:l),vb_help(-l:l) )
c
c-----------------------------------------------------------------------
c
c                __      mm'    *lm'                    __      mm'    *lm' . .
c A_help(G',m) = >     V       A     and B_help(G',m) = >     V       B    <u|u>
c                --m'   l,n,s   G'                      --m'   l,n,s   G'
c
c------------------------------------------------------------------------
      ll = l * (l+1)
      DO m = -l,l
 
        DO mp = -l,l
          va_help(mp) = real( vs_mmp(m,mp,n_l,jsp) ) * invsfct
          vb_help(mp) = va_help(mp) * ddn(l,n,jsp)
        ENDDO

        DO ig = 1, nvv
          ar_help(ig,m) = 0.0
          ai_help(ig,m) = 0.0
          br_help(ig,m) = 0.0
          bi_help(ig,m) = 0.0

          DO mp = -l,l
            ar_help(ig,m) = ar_help(ig,m) + va_help(mp) * 
     +                                      ar(ig,ll+mp,jintsp)
            ai_help(ig,m) = ai_help(ig,m) - va_help(mp) * 
     +                                      ai(ig,ll+mp,jintsp)
            br_help(ig,m) = br_help(ig,m) + vb_help(mp) *
     +                                      br(ig,ll+mp,jintsp)
            bi_help(ig,m) = bi_help(ig,m) - vb_help(mp) *  
     +                                      bi(ig,ll+mp,jintsp)
          ENDDO

        ENDDO

        DO mp = -l,l
          va_help(mp) = aimag( vs_mmp(m,mp,n_l,jsp) ) * invsfct
          vb_help(mp) = va_help(mp) * ddn(l,n,jsp)
        ENDDO
        DO ig = 1, nvv
          DO mp = -l,l
            ar_help(ig,m) = ar_help(ig,m) + va_help(mp) *
     +                                      ai(ig,ll+mp,jintsp)
            ai_help(ig,m) = ai_help(ig,m) + va_help(mp) *
     +                                      ar(ig,ll+mp,jintsp)
            br_help(ig,m) = br_help(ig,m) + vb_help(mp) *
     +                                      bi(ig,ll+mp,jintsp)
            bi_help(ig,m) = bi_help(ig,m) + vb_help(mp) *
     +                                      br(ig,ll+mp,jintsp)
          ENDDO
        ENDDO

      ENDDO
c
      IF (l_ss) nvv = nv(iintsp)
c
c--------------------------------------------
c
c  l,n,s    --        m   lm          m   lm
c H      =  >   A_help   A    + B_help   B
c  G G'     --m       G'  G           G'  G 
c
c--------------------------------------------
      nc = 0
      DO ig = n_rank+1, nvv, n_size
        nc = nc + 1
c
        IF (l_ss) THEN
           IF (n_size > 1) STOP 'U+SS+EV-PAR' 
           IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN     !---> spin-up spin-up part
              igp_m = ig
              chihlp = chi11
              ii = (ig-1)*(ig)/2
           ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN !---> spin-down spin-down part
              igp_m = ig
              chihlp = chi22
              ii = (nv(1)+nlotot+ig-1)*(nv(1)+nlotot+ig)/2 + 
     +              nv(1)+nlotot
           ELSE                                          !---> spin-down spin-up part
              igp_m = nv(1)
              chihlp = chi21
              ii = (nv(1)+nlotot+ig-1)*(nv(1)+nlotot+ig)/2
           ENDIF
        ELSE                                             ! no spin-spiral
           chihlp = cmplx(1.0,0.0)
           igp_m = ig
!          ii = ig * (ig - 1) / 2
           ii = nc*(nc-1)/2*n_size - (nc-1)*(n_size-n_rank-1)
        ENDIF

        DO m = -l,l
          llm = ll + m

          DO igp = 1, igp_m
#ifdef CPP_INVERSION
            aa(ii+igp) = aa(ii+igp) + 
     +                             ar_help(igp,m) * ar(ig,llm,iintsp) +
     +                             br_help(igp,m) * br(ig,llm,iintsp) -
     +                             ai_help(igp,m) * ai(ig,llm,iintsp) -
     +                             bi_help(igp,m) * bi(ig,llm,iintsp) 
#else
            aa(ii+igp) = aa(ii+igp) + chihlp * cmplx(
     +                             ar_help(igp,m) * ar(ig,llm,iintsp) +
     +                             br_help(igp,m) * br(ig,llm,iintsp) -
     +                             ai_help(igp,m) * ai(ig,llm,iintsp) -
     +                             bi_help(igp,m) * bi(ig,llm,iintsp) ,
     +                             ar_help(igp,m) * ai(ig,llm,iintsp) +
     +                             br_help(igp,m) * bi(ig,llm,iintsp) +
     +                             ai_help(igp,m) * ar(ig,llm,iintsp) +
     +                             bi_help(igp,m) * br(ig,llm,iintsp) )
#endif
          ENDDO ! igp 
        ENDDO   ! m
      ENDDO     ! ig = 1, nvv
c--------------------------------------------------------------------------
c
c     local orbitals (if there are any with same l on this atom) 
c                    for the meaning of iilo etc. see slomat
c
c--------------------------------------------------------------------------
      nkvecprevlo = 0
c
      DO lo = 1, nlo(n)
         IF (l.EQ.llo(lo,n)) THEN
         IF (n_size > 1) STOP 'U+LO+EV-PAR'
c
c--->    calculate the contribution to H (G_lo,G): help arrays
c
         ALLOCATE ( cr_help(nvv,-l:l),ci_help(nvv,-l:l) )
         uulo = invsfct * uulon(lo,n)
         dulo = invsfct * dulon(lo,n)
         matel0 = iilo
         locol0 = locol
         ll = l * (l+1)
         DO m = -l,l
           DO mp = -l,l
             va_help(mp) = real( vs_mmp(m,mp,n_l,jsp) ) * uulo
             vb_help(mp) = real( vs_mmp(m,mp,n_l,jsp) ) * dulo
           ENDDO
           DO ig = 1, nvv
             cr_help(ig,m) = 0.0
             ci_help(ig,m) = 0.0
             DO mp = -l,l
              cr_help(ig,m) = cr_help(ig,m) + va_help(mp) *
     +                                      ar(ig,ll+mp,jintsp)
              ci_help(ig,m) = ci_help(ig,m) - va_help(mp) *
     +                                      ai(ig,ll+mp,jintsp)
              cr_help(ig,m) = cr_help(ig,m) + vb_help(mp) *
     +                                      br(ig,ll+mp,jintsp)
              ci_help(ig,m) = ci_help(ig,m) - vb_help(mp) *
     +                                      bi(ig,ll+mp,jintsp)
             ENDDO
           ENDDO
           DO mp = -l,l
             va_help(mp) = aimag( vs_mmp(m,mp,n_l,jsp) ) * uulo
             vb_help(mp) = aimag( vs_mmp(m,mp,n_l,jsp) ) * dulo
           ENDDO
           DO ig = 1, nvv
             DO mp = -l,l
              cr_help(ig,m) = cr_help(ig,m) + va_help(mp) *
     +                                      ai(ig,ll+mp,jintsp)
              ci_help(ig,m) = ci_help(ig,m) + va_help(mp) *
     +                                      ar(ig,ll+mp,jintsp)
              cr_help(ig,m) = cr_help(ig,m) + vb_help(mp) *
     +                                      bi(ig,ll+mp,jintsp)
              ci_help(ig,m) = ci_help(ig,m) + vb_help(mp) *
     +                                      br(ig,ll+mp,jintsp)
             ENDDO
           ENDDO

         ENDDO
c
c--->    now add to H (G_lo,G)
c
         DO m = -l,l
            iilo = matel0
            locol = locol0
            DO nkvec = 1,i_invsf* (2*l+1)
               locol = locol + 1
               IF (mod(locol-1,n_size).EQ.n_rank) THEN
c-t3e
                  DO ig = 1,nvv
                     iilo = iilo + 1
#ifdef CPP_INVERSION
                     aa(iilo) = aa(iilo) + 
     +                    ar_help(ig,m) * real(alo(m,nkvec,lo,iintsp)) +
     +                    br_help(ig,m) * real(blo(m,nkvec,lo,iintsp)) +
     +                    cr_help(ig,m) * real(clo(m,nkvec,lo,iintsp)) -
     +                    ai_help(ig,m) *aimag(alo(m,nkvec,lo,iintsp)) -
     +                    bi_help(ig,m) *aimag(blo(m,nkvec,lo,iintsp)) -
     +                    ci_help(ig,m) *aimag(clo(m,nkvec,lo,iintsp))
#else
                     aa(iilo) = aa(iilo) + cmplx(
     +                    ar_help(ig,m) * real(alo(m,nkvec,lo,iintsp)) +
     +                    br_help(ig,m) * real(blo(m,nkvec,lo,iintsp)) +
     +                    cr_help(ig,m) * real(clo(m,nkvec,lo,iintsp)) -
     +                    ai_help(ig,m) *aimag(alo(m,nkvec,lo,iintsp)) -
     +                    bi_help(ig,m) *aimag(blo(m,nkvec,lo,iintsp)) -
     +                    ci_help(ig,m) *aimag(clo(m,nkvec,lo,iintsp)) ,
     +                    ai_help(ig,m) * real(alo(m,nkvec,lo,iintsp)) +
     +                    bi_help(ig,m) * real(blo(m,nkvec,lo,iintsp)) +
     +                    ci_help(ig,m) * real(clo(m,nkvec,lo,iintsp)) +
     +                    ar_help(ig,m) *aimag(alo(m,nkvec,lo,iintsp)) +
     +                    br_help(ig,m) *aimag(blo(m,nkvec,lo,iintsp)) +
     +                    cr_help(ig,m) *aimag(clo(m,nkvec,lo,iintsp)) )
#endif
                  ENDDO                     ! jump to end of row:
                  iilo = iilo + nkvecprevat + nkvecprevlo + nkvec 
               ENDIF
            ENDDO
         ENDDO
         DEALLOCATE ( cr_help,ci_help )
c
c--->    calculate the contribution to H (G_lo,G_lo)
c
         ALLOCATE ( ax(-l:l,2*(2*llod+1)),wc_help(-l:l),
     +              bx(-l:l,2*(2*llod+1)),wd_help(-l:l),
     +              cx(-l:l,2*(2*llod+1)),we_help(-l:l),
     +              wa_help(-l:l),wb_help(-l:l))

         iilo = matel0
         locol = locol0
         DO nkvec = 1,i_invsf* (2*l+1)
            locol = locol + 1
            IF (mod(locol-1,n_size).EQ.n_rank) THEN
               iilo = iilo + nvv + nkvecprevat
c
c---> with other LOs
c
               DO lop = 1, lo
                  lp = llo(lop,n)
                  IF (l.EQ.lp) THEN
                     dd    = invsfct * ddn(l,n,jsp)
                     ulolo = invsfct * uloulopn(lo,lop,n)
                     nkend = invsfct * (2*lp+1)
                     IF (lop.EQ.lo) nkend = nkvec
                     DO m = -l,l
                       DO mp = -l,l
                         wa_help(mp) = vs_mmp(m,mp,n_l,jsp) * invsfct
                         wb_help(mp) = vs_mmp(m,mp,n_l,jsp) * dd
                         wc_help(mp) = vs_mmp(m,mp,n_l,jsp) * uulo
                         wd_help(mp) = vs_mmp(m,mp,n_l,jsp) * dulo
                         we_help(mp) = vs_mmp(m,mp,n_l,jsp) * ulolo
                       ENDDO
                       DO nkvecp = 1, nkend
                         ax(m,nkvecp) = cmplx(0.0,0.0)
                         bx(m,nkvecp) = cmplx(0.0,0.0)
                         cx(m,nkvecp) = cmplx(0.0,0.0)
                         DO mp = -l,l
                           ax(m,nkvecp) = ax(m,nkvecp) +
     +                  conjg(alo(mp,nkvecp,lop,jintsp)) * wa_help(mp) +
     +                  conjg(clo(mp,nkvecp,lop,jintsp)) * wc_help(mp)
                           bx(m,nkvecp) = bx(m,nkvecp) +
     +                  conjg(blo(mp,nkvecp,lop,jintsp)) * wb_help(mp) +
     +                  conjg(clo(mp,nkvecp,lop,jintsp)) * wd_help(mp)
                           cx(m,nkvecp) = cx(m,nkvecp) +
     +                  conjg(alo(mp,nkvecp,lop,jintsp)) * wc_help(mp) +
     +                  conjg(blo(mp,nkvecp,lop,jintsp)) * wd_help(mp) +
     +                  conjg(clo(mp,nkvecp,lop,jintsp)) * we_help(mp)
                         ENDDO
                       ENDDO
                     ENDDO
c
c---> add to H, special bound for contribution of lo with with itself
c
                     DO nkvecp = 1,nkend  
                        iilo = iilo + 1
                        DO m = -l,l
#ifdef CPP_INVERSION
                           aa(iilo) = aa(iilo) + 
     +                real(alo(m,nkvec,lo,iintsp))* real(ax(m,nkvecp)) +
     +                real(blo(m,nkvec,lo,iintsp))* real(bx(m,nkvecp)) +
     +                real(clo(m,nkvec,lo,iintsp))* real(cx(m,nkvecp)) -
     +               aimag(alo(m,nkvec,lo,iintsp))*aimag(ax(m,nkvecp)) -
     +               aimag(blo(m,nkvec,lo,iintsp))*aimag(bx(m,nkvecp)) -
     +               aimag(clo(m,nkvec,lo,iintsp))*aimag(cx(m,nkvecp)) 
#else
                           aa(iilo) = aa(iilo) + 
     +                           alo(m,nkvec,lo,iintsp) * ax(m,nkvecp) +
     +                           blo(m,nkvec,lo,iintsp) * bx(m,nkvecp) +
     +                           clo(m,nkvec,lo,iintsp) * cx(m,nkvecp) 
#endif
                        ENDDO
                     ENDDO
                  ELSE
                     iilo = iilo + invsfct* (2*lp+1)
                  ENDIF                               ! l = lp
               ENDDO
c
            ENDIF ! T3e
c-t3e
         END DO
         DEALLOCATE ( ax,bx,cx,wa_help,wb_help,wc_help,wd_help,we_help)
         nkvecprevlo = nkvecprevlo + invsfct* (2*l+1)

         ELSE
           lp = invsfct* (2*llo(lo,n)+1)
           iilo = iilo + (nvv + nkvecprevlo + nkvecprevat) * lp
           iilo = iilo + (lp+1) * lp / 2
           nkvecprevlo = nkvecprevlo + lp
           locol = locol + lp
         ENDIF  ! end of if(l=llo(lo,n))
      ENDDO     ! end of lo = 1,nlo loop
      nkvecprevat = nkvecprevat + nkvecprevlo

      DEALLOCATE ( ar_help,ai_help,br_help,bi_help,va_help,vb_help )
c
      END SUBROUTINE u_ham
      END MODULE m_uham
