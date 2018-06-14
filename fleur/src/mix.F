      MODULE m_mix
*************************************************************************
****  mixing of charge densities:                                    ****
****    imix = 0 : linear mixing                                     ****
****    IMIX = 3 : BROYDEN'S FIRST METHOD                            ****
****    IMIX = 5 : BROYDEN'S SECOND METHOD                           ****
****    IMIX = 7 : GENERALIZED ANDERSEN METHOD                       ****
****  implementation to flapw7 ..... R.Pentcheva, D.Vogtenhuber      ****
*************************************************************************
      CONTAINS
      SUBROUTINE mix(
     >               k1d,k2d,k3d,ntypd,ntypsd,jmtd,jspd,nmzd,nmzxyd,n3d,
     >              nrhomfile,film,invs2,invs,nvac,jspins,nq2,nq3,nmzxy,
     >               nmz,z1,delz,ntype,rmt,zatom,dx,rmsh,neq,ntypsy,jri,
     >               nlh,imix,maxiter,alpha,spinf,it,namat,nlhd,n2d,
     >               natd,nop,kv3,nstr,ig,ig2,bmat,sk3,mrot,tau,invtab,
     >               symor,volint,taual,volmts,l_noco,alph,beta,
     >               kimax,igfft,pgfft,ufft,nstr2,area,omtil,sigma,
     >               vol,odi)
c
#include"cpp_double.h"
      USE m_loddop
      USE m_wrtdop
      USE m_brysh1
      USE m_stmix
      USE m_broyden
      USE m_brysh2
      USE m_metric
      USE m_qfix
      USE m_od_types, ONLY : od_inp

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nrhomfile
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,ntypd,ntypsd,jmtd,jspd
      INTEGER, INTENT (IN) :: nmzd,nmzxyd,n3d,nlhd,n2d,nop,natd
      INTEGER, INTENT (IN) :: imix,it,kimax,maxiter,ntype
      INTEGER, INTENT (IN) :: nvac,jspins,nq2,nq3,nmzxy,nmz
      REAL,    INTENT (IN) :: area,alpha,spinf,z1,delz,volint
      REAL,    INTENT (IN) :: sigma,omtil
      LOGICAL, INTENT (IN) :: film,invs,invs2,symor,l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: kv3(3,n3d),nstr(n3d),invtab(nop)
      INTEGER, INTENT (IN) :: nstr2(n2d),nlh(ntypsd),ntypsy(natd)
      INTEGER, INTENT (IN) :: jri(ntypd),neq(ntypd),mrot(3,3,nop)
      REAL,    INTENT (IN) :: dx(ntypd),rmt(ntypd),taual(3,natd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),zatom(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),sk3(n3d),volmts(ntypd)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd)
      CHARACTER*2, INTENT (IN) :: namat(0:103)
c-odim
      REAL,    INTENT (IN) :: vol
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Local Scalars ..
      REAL fix,intfac,vacfac
      INTEGER i,iter,imap,js,mit,nt,nt1,irecl
      INTEGER mmap,mmaph,nmaph,nmap,mapmt,mapvac,mapvac2
      INTEGER iq2,iq3,ivac,imz,n_u,iofl
      LOGICAL lexist,l_ldaU
      CHARACTER*8 dop,iop,name(10)
      INTEGER d1,d10,asciioffset
      CHARACTER (len=5) cdnfile

C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      REAL, ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
      REAL dist(6)
      REAL, ALLOCATABLE :: sm(:),fsm(:) 
c---> off-diagonal part of the density matrix
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      COMPLEX, ALLOCATABLE :: n_mmp(:,:,:,:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC char,sqrt
C     .. External Functions ..
      REAL CPP_BLAS_sdot
      EXTERNAL CPP_BLAS_sdot

c YM: I have exported 'vol' from outside, be aware

c
c     IF (film) THEN
c        vol = 2.0 * z1 * area
c     ELSE
c        vol = omtil
c     ENDIF
c---> In systems without inversions symmetry the interstitial star-
c---> coefficients are complex. Thus twice as many numbers have to be
c---> stored.
      IF (invs) THEN
         intfac = 1.0
      ELSE
         intfac = 2.0
      ENDIF
c---> The corresponding is true for the coeff. of the warping vacuum
c---> density depending on the two dimensional inversion.
      IF (invs2) THEN
         vacfac = 1.0
      ELSE
         vacfac = 2.0
      ENDIF
      mmaph = intfac*n3d + ntypd*(nlhd+1)*jmtd
     +      + vacfac*nmzxyd*(odi%n2d-1)*nvac + nmzd*nvac
      mmap  = mmaph*jspd
c---> in a non-collinear calculations extra space is needed for the
c---> off-diag. part of the density matrix. these coeff. are generally
c---> complex independ of invs and invs2.
      IF (l_noco) mmap = mmap + 2*n3d
     +                  + 2*nmzxyd*(odi%n2d-1)*nvac + 2*nmzd*nvac

      INQUIRE (file='n_mmp_mat',exist=l_ldaU) 
      IF (l_ldaU) THEN
!
! In an LDA+U caclulation, also the density matrix is included in the
! supervectors (sm,fsm) if no mixing factors are in the n_mmp_mat-file
!
        OPEN (69,file='n_mmp_mat',status='old',form='formatted')
        i = 0 
        DO
          READ (69,*,iostat=iofl)
          IF (iofl < 0) EXIT
          i = i + 1
        ENDDO
        IF ( mod(i,14*jspins) == 1 ) THEN    ! was already mixed in u_mix
          n_u = 0
          ALLOCATE ( n_mmp(-3:-3,-3:-3,1,1,2) )
        ELSEIF ( mod(i,28*jspins)== 0 ) THEN ! mix here 
          n_u = i / (28 * jspins )    ! n_u atoms have lda+u applied
          mmap = mmap + 7 * i / 2     ! add 7*7 complex numbers per n_u
          REWIND (69)
          ALLOCATE ( n_mmp(-3:3,-3:3,n_u,jspins,2) )
          READ (69,'(7f20.13)') n_mmp(:,:,:,:,:)
        ELSE
          STOP 'strange n_mmp_mat-file...'
        ENDIF
        CLOSE (69)
      ELSE
        n_u=0
        ALLOCATE ( n_mmp(-3:-3,-3:-3,1,1,2) )
      ENDIF
!
      ALLOCATE (sm(mmap),fsm(mmap))

      ALLOCATE (qpw(n3d,jspd),rhtxy(nmzxyd,odi%n2d-1,2,jspd),
     &          rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd) )

      IF (l_noco) THEN
        ALLOCATE (cdom(n3d),cdomvz(nmzd,2),
     +                      cdomvxy(nmzxyd,odi%n2d-1,2))
      ELSE
        ALLOCATE ( cdom(1),cdomvz(1,1),cdomvxy(1,1,1) )
      ENDIF
c
      INQUIRE (file='broyd.'//char(imix+48),exist=lexist)
      DO i = 1,6
            dist(i) = 0.0
      END DO
      mit = 0
c---> detremine type of mixing:
c---> imix=0:straight, imix=o broyden first, imix=5:broyden second
c---> imix=:generalozed anderson mixing
      IF (imix.EQ.0) THEN
         WRITE (16,FMT='(a,2f10.5)') 'STRAIGHT MIXING',alpha
      ELSEIF (imix.EQ.3) THEN
         IF ( .not.lexist) mit = 1
         WRITE (16,FMT='(a,f10.5)') 'BROYDEN FIRST MIXING',alpha
      ELSEIF (imix.EQ.5) THEN
         IF (.not.lexist) mit = 1
         WRITE (16,FMT='(a,f10.5)') 'BROYDEN SECOND MIXING',alpha
      ELSEIF (imix.EQ.7) THEN
         IF (.not.lexist) mit = 1
         WRITE (16,FMT='(a,f10.5)') 'ANDERSON GENERALIZED',alpha
      ELSE
         STOP 'mix: imix =/= 0,3,5,7 '
      END IF
c
      IF (jspins.EQ.2.and.imix.ne.0) THEN
         WRITE(6,'(''WARNING : for QUASI-NEWTON METHODS SPINF=1'')')
      END IF

c--->    generate name of file to hold the results of this iteration
!gs+  changed cdn filename gen. to cdn//tens//digits
!gs   fixed cdn filehandle to 71, otherwise conflict after it>=18
!gs   plus: you don't want to have open so many fh's
!gs    construct 2char rep of iter --> '01'-->'Z9'
      d1 = mod(it,10)
      d10 = int( (it + 0.5)/10 )
      asciioffset = iachar('1')-1
      IF ( d10.GE.10 ) asciioffset = iachar('7')
      cdnfile = 'cdn'//achar(d10+asciioffset)//achar(d1+iachar('1')-1)
!      WRITE (*,*) d1,d10,'>',cdnfile,'<' !gsdbg
      IF (.NOT.l_noco) THEN
         nt1 = 72
         OPEN (nt1,file=cdnfile,form='unformatted',status='unknown')
         REWIND nt1
      ENDIF
!gs-

c---> reload densities of current iteration
      IF (l_noco) THEN
c---> initialize arrays for the off-diagonal part of the density matrix
         cdom(:) = cmplx(0.0,0.0)
         IF (film) THEN
            cdomvz(:,:) = cmplx(0.0,0.0)
            cdomvxy(:,:,:) = cmplx(0.0,0.0)
         END IF

c--->    reload input density matrix from file rhomat_inp
         OPEN (nrhomfile,FILE='rhomat_inp',FORM='unformatted',
     +         STATUS='unknown')
c--->    first the diagonal elements of the density matrix
         CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nrhomfile,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c--->    and then the off-diagonal part
         READ (nrhomfile,END=150,ERR=50) (cdom(iq3),iq3=1,nq3)
         IF (film) THEN
            READ (nrhomfile,END=75,ERR=50) ((cdomvz(imz,ivac),
     +                                  imz=1,nmz),ivac=1,nvac)
            READ (nrhomfile,END=75,ERR=50) (((cdomvxy(imz,iq2-1,ivac),
     +                        imz=1,nmzxy),iq2=2,odi%nq2),ivac=1,nvac)
         ENDIF
         GOTO 150
 50      WRITE(6,*)'rhodirgen: ERROR: Problems while reading density'
         WRITE(6,*)'matrix from file rhomat_inp.'
         STOP 'rhomatdir: ERROR while reading file rhomat_inp'
 75      WRITE(6,*)'rhomatdir: ERROR: reached end of file rhomat_inp'
         WRITE(6,*)'while reading the vacuum part of the off-diagonal'
         WRITE(6,*)'element of the desity matrix.'
         STOP 'rhomatdir: ERROR while reading file rhomat_inp'
 150     CLOSE (nrhomfile)
      ELSE
         nt = 71                 !gs see above
         OPEN (nt,file='cdn1',form='unformatted',status='old')
         REWIND nt
         CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c--->      write density to file for storage
         CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt1,
     >            iop,dop,iter,rho,qpw,rht,rhtxy,name)
      ENDIF
c
c--->  put input charge density into arrays sm 
c      in the spin polarized case the arrays consist of 
c      spin up and spin down densities

      CALL brysh1(
     >            jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap,n2d,natd,
     >            l_noco,nq3,ntype,nlh,ntypsy,jri,nmzxy,nmz,
     >            film,invs,invs2,nvac,nq2,neq,intfac,vacfac,
     >            qpw,rho,rht,rhtxy,cdom,cdomvz,cdomvxy,
     >            n_u,n_mmp(-3,-3,1,1,1),odi,
     <            nmap,nmaph,mapmt,mapvac,mapvac2,sm) 
c
c     load output charge density
c
      IF (l_noco) THEN
c--->    reload output density matrix from file rhomat_out
         OPEN (nrhomfile,FILE='rhomat_out',FORM='unformatted',
     +         STATUS='unknown')
c--->    first the diagonal elements of the density matrix
         CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nrhomfile,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c--->    and then the off-diagonal part
         READ (nrhomfile) (cdom(iq3),iq3=1,nq3)
         IF (film) THEN
            READ (nrhomfile) ((cdomvz(imz,ivac),imz=1,nmz),ivac=1,nvac)
            READ (nrhomfile) (((cdomvxy(imz,iq2-1,ivac),imz=1,nmzxy),
     +                                       iq2=2,odi%nq2),ivac=1,nvac)
         ENDIF
         CLOSE (nrhomfile)
      ELSE
         CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c--->      write density to file for storage
         CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt1,
     >            iop,dop,iter,rho,qpw,rht,rhtxy,name)
         CLOSE(nt1)
      ENDIF
c
c--->  put output charge density into arrays fsm 
c
      CALL brysh1(
     >            jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap,n2d,natd,
     >            l_noco,nq3,ntype,nlh,ntypsy,jri,nmzxy,nmz,
     >            film,invs,invs2,nvac,nq2,neq,intfac,vacfac,
     >            qpw,rho,rht,rhtxy,cdom,cdomvz,cdomvxy,
     >            n_u,n_mmp(-3,-3,1,1,2),odi,
     <            nmap,nmaph,mapmt,mapvac,mapvac2,fsm) 
c
c --> store fsm - sm the difference on fsm
c      
      DO imap = 1,nmap
         fsm(imap) = fsm(imap) - sm(imap)
      END DO
c
c open files for broyden
c
      irecl=(2*nmap+1)*8
      IF (imix.GE.3) THEN
         OPEN (57,file='broyd',form='unformatted',status='unknown')
         OPEN (59,file='broyd.'//char(imix+48),access='direct',
     +             recl=irecl,form='unformatted',status='unknown')
      END IF
c
c----->  mixing of the densities
c
      IF (imix.eq.0) THEN
         CALL stmix(
     >              mmap,jspd,n_u,
     >              nmap,nmaph,jspins,l_noco,alpha,spinf,fsm,
     =              sm)
      ELSE
         CALL broyden(
     >                z1,jspd,n3d,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,
     >                nmaph,mapmt,mapvac2,n2d,k1d,k2d,k3d,natd,
     >                nmap,jspins,l_noco,imix,maxiter,alpha,fsm,ufft,
     >                kimax,igfft,pgfft,nq3,omtil,odi,
     >                ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >                film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,
     =                mit,sm)
      END IF
c     call timestamp(-3)
c----->  load output densities
c
      CALL brysh2(
     >            jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap,n2d,natd,
     >            l_noco,nq3,ntype,nlh,ntypsy,jri,nmzxy,nmz,
     >            film,invs,invs2,nvac,nq2,neq,nmap,sm,n_u,
     <            n_mmp(-3,-3,1,1,1),odi,
     <            qpw,rho,rht,rhtxy,cdom,cdomvz,cdomvxy) 

c
c----->  calculate the distance of charge densities
c
c     induce metric in fsm use sm as an output array:
c     |sm> = w |fsm>
c
      CALL metric(
     >     z1,jspd,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,nmaph,
     >     mapmt,mapvac2,jspins,l_noco,fsm,ufft,kimax,igfft,
     >     pgfft,nq3,omtil,k1d,k2d,k3d,n2d,n3d,natd,
     >     ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >     film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,odi,
     <     sm) 
c
c     calculate the charge density distance for each spin
c
      DO 30 js = 1,jspins
         dist(js) = CPP_BLAS_sdot(nmaph,fsm(nmaph*(js-1)+1),1,
     +                        sm(nmaph*(js-1)+1),1)
         WRITE (16,FMT=7900) js,iter,1000*sqrt(ABS(dist(js)/vol))
         WRITE (6,FMT=7900) js,iter,1000*sqrt(ABS(dist(js)/vol))
 30   CONTINUE
      IF (l_noco) dist(6) = CPP_BLAS_sdot((nmap-2*nmaph),
     +                             fsm(nmaph*2+1),1,sm(nmaph*2+1),1)
      IF (l_noco) WRITE (6,FMT=7900) 3,iter,1000*sqrt(ABS(dist(6)/vol))
c
c     calculate the distance of total charge and spin density
c     |rho/m(o) - rho/m(i)| = |rh1(o) -rh1(i)|+ |rh2(o) -rh2(i)| +/_
c                             +/_2<rh2(o) -rh2(i)|rh1(o) -rh1(i)>
c
      IF (jspins.EQ.2) THEN
         dist(3) = CPP_BLAS_sdot(nmaph,fsm,1,sm(nmaph+1),1)
         dist(4) = dist(1) + dist(2) + 2.0e0*dist(3)
         dist(5) = dist(1) + dist(2) - 2.0e0*dist(3)
         WRITE (16,FMT=8000) iter,1000*sqrt(ABS(dist(4)/vol))
         WRITE (16,FMT=8010) iter,1000*sqrt(ABS(dist(5)/vol))
         WRITE (6,FMT=8000) iter,1000*sqrt(ABS(dist(4)/vol))
         WRITE (6,FMT=8010) iter,1000*sqrt(ABS(dist(5)/vol))
         ! dist/vol should always be >= 0 ,
         ! but for dist=0 numerically you might obtain dist/vol < 0
         ! (e.g. when calculating non-magnetic systems with jspins=2).
      END IF
      DEALLOCATE (sm,fsm)
c
c----> output of mixed densities
c
c     ---> fix the new density
      CALL qfix(
     >          k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >          nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,nmzxy,n2d,
     >          ntype,neq,volmts,taual,z1,vol,volint,nq2,invtab,
     >          symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >          nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     >          qpw,rhtxy,rho,rht,odi,
     <          fix)

      iter = iter + 1
      IF (l_noco) THEN
c
c--->    write mixed density to file rhomat_out
c
         OPEN (nrhomfile,FILE='rhomat_inp',FORM='unformatted',
     +         STATUS='unknown')
c--->    first the diagonal elements of the density matrix
         CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nrhomfile,
     >            'input   ','density ',iter,rho,qpw,rht,rhtxy,name)
c--->    and then the off-diagonal part
         WRITE (nrhomfile) (cdom(iq3),iq3=1,nq3)
         IF (film) THEN
            WRITE (nrhomfile) ((cdomvz(imz,ivac),imz=1,nmz),ivac=1,nvac)
            WRITE (nrhomfile) (((cdomvxy(imz,iq2-1,ivac),imz=1,nmzxy),
     +                                      iq2=2,odi%nq2),ivac=1,nvac)
         ENDIF
         CLOSE (nrhomfile)
      ELSE
c
c--->    write new density onto unit 71 (overwrite)
c
         REWIND nt
         CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt,
     >            'input   ','density ',iter,rho,qpw,rht,rhtxy,name)
         close (nt)
      ENDIF
      DEALLOCATE ( cdom,cdomvz,cdomvxy )
      IF ( n_u > 0 ) THEN
        OPEN (69,file='n_mmp_mat',status='replace',form='formatted')
        WRITE (69,'(7f20.13)') n_mmp(:,:,:,:,1)
        CLOSE (69)
      ENDIF
      DEALLOCATE (n_mmp)

      IF (imix.GT.0) THEN
         CLOSE (57)
         CLOSE (59)
      END IF
 7900 FORMAT (/,'----> distance of charge densities for spin ',i2,' 
     +       it=',i3,':',f11.6,' me/bohr**3')
 8000 FORMAT (/,'----> distance of charge densities for it=',i3,':',
     +       f11.6,' me/bohr**3')
 8010 FORMAT (/,'----> distance of spin densities for it=',i3,':',
     +         f11.6,' me/bohr**3')
 8020 FORMAT (4d25.14)
 8030 FORMAT (10i10)

      DEALLOCATE (qpw,rhtxy,rho,rht)

      END SUBROUTINE mix
      END MODULE m_mix
