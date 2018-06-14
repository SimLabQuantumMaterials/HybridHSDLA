      MODULE m_tlmplm
c*********************************************************************
c     sets up the t(l'm',lm) matrices for each atom type. these
c     matrices are the k-point independent quantities needed to
c     determine the hamiltonian contributions coming from
c     the spheres.  the results are written to direct access
c     unit28 labelled by atom type number (this allows for
c     multitasking).  only non=zero elements are calculated
c     and an index array is also set up.
c               m. weinert     1986
!---------------------------------------------------------------------
!     introduced lh0 and vr0 for gw-calculations (see below)
c*********************************************************************
      CONTAINS
      SUBROUTINE tlmplm(
     >                  memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd,
     >                  ntype,dx,rmsh,jri,lmax,ntypsy,nwdd,natd,
     >                  lnonsph,lmd,lmplmd,clnu,mlh,nmem,llh,nlh,neq,
     >                  jsp,nw,nb0,irank,mlotot,mlolotot,
     >                  vr,vr0,gwc,lh0,el,secvar,
     >                  nlod,llod,loplod,ello,llo,nlo,lo1l,l_dulo,
     >                  ulo_der,l_f,lmaxb,n_u,lda_u,vs_mmp,
     <                  tuu,tud,tdu,tdd,ind,tuulo,tdulo,tuloulo,
     <                  us,dus,uds,duds,ddn,
     <                  ulos,dulos,uulon,dulon,uloulopn)

      USE m_intgr, ONLY : intgr3
      USE m_radflo
      USE m_radfun
      USE m_tlo
      USE m_gaunt, ONLY: gaunt1
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd
      INTEGER, INTENT (IN) :: lmd,lmplmd,jsp,nb0,nw,ntype,nwdd,irank
      INTEGER, INTENT (IN) :: nlod,llod,loplod,natd,mlotot,mlolotot
      INTEGER, INTENT (IN) :: lmaxb,n_u,gwc,lh0 ! normally lh0 = 1; for second round (gwc=2) of gw-calcs lh0 = 0
      LOGICAL, INTENT (IN) :: secvar,l_f
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd),lda_u(ntypd)
      INTEGER, INTENT (IN) :: lo1l(0:llod,ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),lmax(ntypd),ntypsy(natd)
      INTEGER, INTENT (IN) :: lnonsph(ntypd),ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd,ntypd)   ! this is for the integrals 
      REAL,    INTENT (IN) :: vr0(jmtd,ntypd)         ! just for radial functions
      REAL,    INTENT (IN) :: el(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: uds(0:lmaxd,ntypd),us(0:lmaxd,ntypd)
      REAL,    INTENT (OUT):: duds(0:lmaxd,ntypd),dus(0:lmaxd,ntypd)
      REAL,    INTENT (OUT):: ulos(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (OUT):: uulon(nlod,ntypd),dulon(nlod,ntypd)
      REAL,    INTENT (OUT):: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (OUT):: ddn(0:lmaxd,ntypd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      COMPLEX, INTENT (OUT):: tdd(0:lmplmd,ntypd),tdu(0:lmplmd,ntypd)
      COMPLEX, INTENT (OUT):: tud(0:lmplmd,ntypd),tuu(0:lmplmd,ntypd)
      COMPLEX, INTENT (OUT):: tdulo(0:lmd,-llod:llod,mlotot)
      COMPLEX, INTENT (OUT):: tuulo(0:lmd,-llod:llod,mlotot)
      COMPLEX, INTENT (OUT):: tuloulo(-llod:llod,-llod:llod,mlolotot)
      INTEGER, INTENT (OUT):: ind(0:lmd,0:lmd,ntypd)
C     ..
C     .. Local Scalars ..
      COMPLEX cil,ci
      REAL temp,wronk
      INTEGER i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl,lmplm,lmx,lmxx,lp
      INTEGER lp1,lpl,m,mem,mems,mp,mu,n,nh,noded,nodeu,nrec,nsym,na,ne
      INTEGER mlo, mlolo, iu
C     ..
C     .. Local Arrays ..
      INTEGER indt(0:lmplmd)
      REAL dvd(0:lmaxd*(lmaxd+3)/2,lh0:nlhd)
      REAL dvu(0:lmaxd*(lmaxd+3)/2,lh0:nlhd)
      REAL uvd(0:lmaxd*(lmaxd+3)/2,lh0:nlhd)
      REAL uvu(0:lmaxd*(lmaxd+3)/2,lh0:nlhd)
      REAL f(jmtd,2,0:lmaxd),g(jmtd,2,0:lmaxd),x(jmtd)
      REAL flo(jmtd,2,nlod)
      REAL uuilon(nlod,ntypd),duilon(nlod,ntypd)
      REAL ulouilopn(nlod,nlod,ntypd)

#if (defined(CPP_MPI) && !defined(CPP_T90))
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,max,mod
C     ..
      ci = cmplx(0.0,1.0)
      na = 1
      mlo = 1 ; mlolo = 1

      tdulo(:,:,:) = cmplx(0.0,0.0)
      tuulo(:,:,:) = cmplx(0.0,0.0)
      tuloulo(:,:,:) = cmplx(0.0,0.0)

      DO 210 n = 1,ntype
c
c--->    generate the wavefunctions for each l
c
         IF (irank.eq.0) WRITE (6,FMT=8000) n
         IF (gwc==2) WRITE (14) rmsh(1:jri(n),n),lmax(n)
         DO 10 l = 0,lmax(n)
            CALL radfun(
     >                  l,el(l,n,jsp,nw),vr0(1,n),jri(n),rmsh(1,n),
     >                  dx(n),jmtd,
     <                  f(1,1,l),g(1,1,l),us(l,n),dus(l,n),uds(l,n),
     <                  duds(l,n),ddn(l,n),nodeu,noded,wronk)
            IF (irank.eq.0) WRITE (6,FMT=8010) l,el(l,n,jsp,nw),us(l,n),
     +            dus(l,n),nodeu,uds(l,n),duds(l,n),noded,ddn(l,n),wronk
            IF (gwc==2) WRITE (14) f(1:jri(n),1,l),g(1:jri(n),1,l),
     &                             f(1:jri(n),2,l),g(1:jri(n),2,l)
   10    CONTINUE
 8000    FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',
     +          /,t32,'radial function',t79,'energy derivative',/,t3,
     +          'l',t8,'energy',t26,'value',t39,'derivative',t53,
     +          'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,
     +          'norm',t119,'wronskian')
 8010    FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)
c
c--->   generate the extra wavefunctions for the local orbitals,
c--->   if there are any.
c
         IF (nlo(n).GE.1) THEN
            CALL radflo(
     >                  ntypd,nlod,jspd,jmtd,lmaxd,n,jsp,ello(1,1,jsp),
     >                  vr0(1,n),jri(n),rmsh(1,n),dx(n),
     >                  f,g,llo,nlo,l_dulo(1,n),irank,ulo_der,
     <                  ulos,dulos,uulon,dulon,uloulopn,
     <                  uuilon,duilon,ulouilopn,flo)
            IF (gwc==2) THEN
              DO i=1,nlo(n)
                WRITE (14) flo(1:jri(n),1,i),
     &                     flo(1:jri(n),2,i)
              ENDDO
            ENDIF
         END IF

           nsym = ntypsy(na)
           nh = nlh(nsym)
c
c--->    generate the irreducible integrals (u(l'):v(lamda,nu:u(l))
c--->    for l' .ge. l, but only thos that will contribute
c
         DO 80 lp = 0,lmax(n)
            lp1 = (lp* (lp+1))/2
            DO 70 l = 0,lp
               lpl = lp1 + l
c--->    loop over non-spherical components of the potential: must
c--->    satisfy the triangular conditions and that l'+l+lamda even
c--->    (conditions from the gaunt coefficient)
               DO 60 lh = lh0, nh
                  lamda = llh(lh,nsym)
                  lmin = lp - l
                  lmx = lp + l
                  IF ((mod(lamda+lmx,2).EQ.1) .OR. (lamda.LT.lmin) .OR.
     +                (lamda.GT.lmx)) THEN
                     uvu(lpl,lh) = 0.0
                     dvd(lpl,lh) = 0.0
                     uvd(lpl,lh) = 0.0
                     dvu(lpl,lh) = 0.0
                  ELSE
                     DO 20 i = 1,jri(n)
                        x(i) = (f(i,1,lp)*f(i,1,l)+f(i,2,lp)*f(i,2,l))*
     +                         vr(i,lh,n)
   20                CONTINUE
                     CALL intgr3(x,rmsh(1,n),dx(n),jri(n),temp)
                     uvu(lpl,lh) = temp
                     DO 30 i = 1,jri(n)
                        x(i) = (g(i,1,lp)*f(i,1,l)+g(i,2,lp)*f(i,2,l))*
     +                         vr(i,lh,n)
   30                CONTINUE
                     CALL intgr3(x,rmsh(1,n),dx(n),jri(n),temp)
                     dvu(lpl,lh) = temp
                     DO 40 i = 1,jri(n)
                        x(i) = (f(i,1,lp)*g(i,1,l)+f(i,2,lp)*g(i,2,l))*
     +                         vr(i,lh,n)
   40                CONTINUE
                     CALL intgr3(x,rmsh(1,n),dx(n),jri(n),temp)
                     uvd(lpl,lh) = temp
                     DO 50 i = 1,jri(n)
                        x(i) = (g(i,1,lp)*g(i,1,l)+g(i,2,lp)*g(i,2,l))*
     +                         vr(i,lh,n)
   50                CONTINUE
                     CALL intgr3(x,rmsh(1,n),dx(n),jri(n),temp)
                     dvd(lpl,lh) = temp
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
c--->    generate the various t(l'm',lm) matrices for l'm'.ge.lm
         tuu(:,n) = cmplx(0.0,0.0)
         tdd(:,n) = cmplx(0.0,0.0)
         tud(:,n) = cmplx(0.0,0.0)
         tdu(:,n) = cmplx(0.0,0.0)
         indt(:) = 0
c--->    loop over l'm'
         DO 140 lp = 0,lmax(n)
            lp1 = (lp* (lp+1))/2
            DO 130 mp = -lp,lp
               lmp = lp* (lp+1) + mp
               lmpl = (lmp* (lmp+1))/2
c--->    loop over lattice harmonics
               DO 120 lh = lh0, nh
                  lamda = llh(lh,nsym)
                  lmin0 = abs(lp-lamda)
                  IF (lmin0.GT.lp) GO TO 120
c-->     ensure l+l'+lamda even
                  lmxx = lp - mod(lamda,2)
                  mems = nmem(lh,nsym)
                  DO 110 mem = 1,mems
                     mu = mlh(mem,lh,nsym)
                     m = mp - mu
                     lmin = max(lmin0,abs(m))
                     l2 = abs(lmxx-lmin)
                     lmin = lmin + mod(l2,2)
                     DO 100 l = lmin,lmxx,2
                        lm = l* (l+1) + m
                        IF (lm.GT.lmp) GO TO 100
                        lpl = lp1 + l
                        lmplm = lmpl + lm
                        cil = ((ci** (l-lp))*clnu(mem,lh,nsym))*
     +                        gaunt1(lp,lamda,l,mp,mu,m,lmaxd)
                        tuu(lmplm,n) = tuu(lmplm,n) + cil*uvu(lpl,lh)
                        tdd(lmplm,n) = tdd(lmplm,n) + cil*dvd(lpl,lh)
                        tud(lmplm,n) = tud(lmplm,n) + cil*uvd(lpl,lh)
                        tdu(lmplm,n) = tdu(lmplm,n) + cil*dvu(lpl,lh)
                        indt(lmplm) = 1
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
c--->    include diagonal terms from muffin-tin hamiltonian
c--->   only, if secvar = .false. If second variation is          
c--->   is switched on, the tlmplm matrix contains only the       
c--->   off-diagonal terms                                       
C-gu deleted: 20.11.97
c--->    set up mapping array
         DO 200 lp = 0,lmax(n)
            DO 190 mp = -lp,lp
               lmp = lp* (lp+1) + mp
               DO 180 l = 0,lmax(n)
                  DO 170 m = -l,l
                     lm = l* (l+1) + m
                     IF (lmp.GE.lm) THEN
                        lmplm = (lmp* (lmp+1))/2 + lm
                        IF (indt(lmplm).NE.0) THEN
                           ind(lmp,lm,n) = lmplm
                        ELSE
                           ind(lmp,lm,n) = -9999
                        END IF
                     ELSE
                        lmplm = (lm* (lm+1))/2 + lmp
                        IF (indt(lmplm).NE.0) THEN
                           ind(lmp,lm,n) = -lmplm
                        ELSE
                           ind(lmp,lm,n) = -9999
                        END IF
                     END IF
  170             CONTINUE
  180          CONTINUE
  190       CONTINUE
  200    CONTINUE
c
c--->   set up the t-matrices for the local orbitals,
c--->   if there are any
         IF (nlo(n).GE.1) THEN
            CALL tlo(
     >               lmaxd,ntypd,jspd,nwdd,jmtd,nlod,llod,nlhd,ntypsd,
     >               lmd,loplod,memd,nlh,llh,nmem,mlh,clnu,ntypsy,
     >               jsp,n,lmax(n),jri(n),rmsh(1,n),dx(n),el,ello,
     >               lh0,secvar,vr(1,0,n),natd,na,l_dulo,ulo_der,
     >               flo,f,g,nlo,llo,lo1l,uulon,dulon,uloulopn,
     >               uuilon,duilon,ulouilopn,
     <               tuulo(0,-llod,mlo),tdulo(0,-llod,mlo),
     <               tuloulo(-llod,-llod,mlolo))
            mlo = mlo + nlo(n)
            mlolo = mlolo + nlo(n)*(nlo(n)+1)/2
         ENDIF
         na = na + neq(n)
  210 ENDDO

      IF (l_f) THEN                  !--->    store on unit 28/38
         mlo = 1 ; mlolo = 1 ; iu = 0
         DO n = 1, ntype
         nrec = nb0 + n
         IF (nlo(n).GE.1) THEN
#if (defined(CPP_MPI) && !defined(CPP_T90))
           IF (irank.eq.0) THEN
#endif
           IF (lda_u(n).GE.0) THEN
             iu = iu + 1
             IF (jsp.EQ.1) WRITE (28,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n),
     +                                tuulo(:,:,mlo:mlo+nlo(n)-1),
     +                                tdulo(:,:,mlo:mlo+nlo(n)-1),
     +             tuloulo(:,:,mlolo:mlolo+nlo(n)*(nlo(n)+1)/2-1),
     +                                vs_mmp(:,:,iu)
             IF (jsp.EQ.2) WRITE (38,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n),
     +                                tuulo(:,:,mlo:mlo+nlo(n)-1),
     +                                tdulo(:,:,mlo:mlo+nlo(n)-1),
     +             tuloulo(:,:,mlolo:mlolo+nlo(n)*(nlo(n)+1)/2-1),
     +                                vs_mmp(:,:,iu)
           ELSE
             IF (jsp.EQ.1) WRITE (28,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n),
     +                                tuulo(:,:,mlo:mlo+nlo(n)-1),
     +                                tdulo(:,:,mlo:mlo+nlo(n)-1),
     +              tuloulo(:,:,mlolo:mlolo+nlo(n)*(nlo(n)+1)/2-1)
             IF (jsp.EQ.2) WRITE (38,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n),
     +                                tuulo(:,:,mlo:mlo+nlo(n)-1),
     +                                tdulo(:,:,mlo:mlo+nlo(n)-1),
     +              tuloulo(:,:,mlolo:mlolo+nlo(n)*(nlo(n)+1)/2-1)
           ENDIF 
#if (defined(CPP_MPI) && !defined(CPP_T90))
           ENDIF
#endif
           mlo = mlo + nlo(n)
           mlolo = mlolo + nlo(n)*(nlo(n)+1)/2
         ELSE
#if (defined(CPP_MPI) && !defined(CPP_T90))
           IF (irank.eq.0) THEN
#endif
           IF (lda_u(n).GE.0) THEN
             iu = iu + 1
             IF (jsp.EQ.1) WRITE (28,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n),
     +                                vs_mmp(:,:,iu)
             IF (jsp.EQ.2) WRITE (38,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n),
     +                                vs_mmp(:,:,iu)
           ELSE
             IF (jsp.EQ.1) WRITE (28,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n)
             IF (jsp.EQ.2) WRITE (38,rec=nrec) tuu(:,n),tdd(:,n),
     +                                tud(:,n),tdu(:,n),ind(:,:,n)
           ENDIF
#if (defined(CPP_MPI) && !defined(CPP_T90))
           ENDIF
#endif
         ENDIF
         ENDDO
      ENDIF
#if (defined(CPP_MPI) && !defined(CPP_T90))
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      END SUBROUTINE tlmplm
      END MODULE m_tlmplm
