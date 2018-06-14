      MODULE m_eigenso
c
c*********************************************************************
c     sets up and solves the spin-orbit eigenvalue problem in the
c     second variation procedure.
c
c     way: takes e.v. and e.f. from previous scalar-rel. calc.
c     makes spin-orbit matrix elements solves e.v. and put it on 'eig'
c
c     Tree:  eigenso-|- loddop
c                    |- spnorb  : sets up s-o parameters 
c                    |    |- soinit - sorad  : radial part
c                    |    |- sgml            : diagonal angular parts
c                    |    |- anglso          : non-diagonal -"-
c                    |
c                    |- alineso : sets up and solves e.v. problem
c                         |- hsohelp
c                         |- hsoham
c
c**********************************************************************
c
      CONTAINS
      SUBROUTINE eigenso(
     > jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,nlod,
     > nwdd,neigd,nvd,lmaxd,nop,natd,ntypsd,llod,
     > ulo_der,l_dulo,
     > invsat,invsatnr,lapw_l,bbmat,l_noco,alph,beta,
     > nbasfcn,film,invs,invs2,nq2,nq3,nvac,jspins,nwd,
     > irecl1,nu,irank,isize,irecl,ntype,soc_opt,nlo,llo,
     > thetaIn,phiIn,omtil,neq,lmax,mrot,ngopr,nkpt,nlh,
     > jri,ntypsy,bmat,rmt,dx,rmsh,taual,invtab,fname,
     > odi,ods,eonly,TYP_REAL)
c
      USE m_od_types, ONLY : od_inp, od_sym
#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : openeig,read_kptenpar,writeeigc,closeeig
#endif
      USE m_spnorb 
      USE m_alineso
      USE m_loddop
      USE m_cputime
      USE m_outtime
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
      USE m_mpi_col_eig
#endif
      IMPLICIT NONE
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
#endif
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd
      INTEGER, INTENT (IN) :: nwdd,neigd,nvd,lmaxd,nop,natd,llod
      INTEGER, INTENT (IN) :: ntype,nu,irecl,nq2,nq3,nvac,jspins
      INTEGER, INTENT (IN) :: irank,isize,ntypsd,irecl1,nlod,nwd 
      INTEGER, INTENT (IN) :: nbasfcn
      REAL,    INTENT (IN) :: thetaIn,phiIn,omtil
      LOGICAL, INTENT (IN) :: film,invs,invs2,l_noco,eonly 
      INTEGER, INTENT (IN) :: TYP_REAL
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),nkpt(nwdd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),jri(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd),lapw_l(ntype)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: bmat(3,3),dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: rmt(ntypd),taual(3,natd),bbmat(3,3)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd)
      LOGICAL, INTENT (IN) :: soc_opt(ntypd+2) 
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      CHARACTER(len=12),INTENT (IN) :: fname(3)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
#ifdef CPP_MPI
      INTEGER ierr(3),stt(MPI_STATUS_SIZE)
#endif
      INTEGER i,j,nk,nk0,nw,jspin,nv,iter,nrec,nk1,n,nlotot,l
      INTEGER n_loc,n_plus,i_plus,n_end,nsz,nmat,nk00,nv1(2)
      REAL    cp_spnorb,cp_aline,t1,t2,cp_help,cp_ham,cp_diag,wtkpt,wk
      REAL    theta, phi, pi
      LOGICAL l_all
      LOGICAL oldcv
C     ..
C     .. Local Arrays..
      CHARACTER*8 dop,iop,name(10)
      CHARACTER*3 chntype

      INTEGER, ALLOCATABLE :: k1(:),k2(:),k3(:),kveclo(:)
      REAL,    ALLOCATABLE :: rk(:),bkpt(:)
      REAL,    ALLOCATABLE :: evac(:,:,:),ello(:,:,:),el(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopdp(:,:,:,:),rsopdpd(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopp(:,:,:,:),rsoppd(:,:,:,:),ulos(:,:,:)
      REAL,    ALLOCATABLE :: us(:,:,:),dus(:,:,:),uds(:,:,:),ddn(:,:,:)
      REAL,    ALLOCATABLE :: duds(:,:,:),eig_so(:),dulos(:,:,:)
      REAL,    ALLOCATABLE :: uulon(:,:,:),dulon(:,:,:),rsoplop(:,:,:,:)
      REAL,    ALLOCATABLE :: rsoplopd(:,:,:,:),rsopdplo(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopplo(:,:,:,:),rsoploplop(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: zso(:,:,:),soangl(:,:,:,:,:,:)

      REAL,    ALLOCATABLE :: vz(:,:,:),vr(:,:,:,:)
      COMPLEX, ALLOCATABLE :: vzxy(:,:,:,:),vpw(:,:)

      INTEGER, ALLOCATABLE :: impi(:)
      REAL,    ALLOCATABLE :: rmpi(:)
C  ..

      oldcv= .true. ! .false.
      ! oldcv=T <=> theta and phi specified according to old convention

c
      IF (oldcv) THEN
        theta= thetaIn
        phi= phiIn
      ELSE
        pi= 4.*ATAN(1.)
        theta= -thetaIn
        phi=   phiIn+pi
        ! now the definition of rotation matrices
        ! is equivalent to the def in the noco-routines
      ENDIF
c
      cp_spnorb = 0.0
      cp_aline  = 0.0
      cp_help   = 0.0
      cp_ham    = 0.0
      cp_diag   = 0.0
c
c load potential from file pottot (=unit 8)
c
      ALLOCATE ( vz(nmzd,2,jspd),vr(jmtd,0:nlhd,ntypd,jspd),
     +           vzxy(nmzxyd,odi%n2d-1,2,jspd),vpw(nq3,jspd) )

      OPEN (nu,file='pottot',form='unformatted',status='old')
      CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nu,natd,neq,
     <            iop,dop,iter,vr,vpw,vz,vzxy,name)
      CLOSE(nu)

      DEALLOCATE ( vz,vzxy,vpw )
c
c     determine the total number of lo's : nlotot
c
      nlotot = 0
      DO n = 1, ntype
        DO l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        ENDDO
      ENDDO
      ALLOCATE(  us(0:lmaxd,ntypd,jspd), dus(0:lmaxd,ntypd,jspd),
     +          uds(0:lmaxd,ntypd,jspd),duds(0:lmaxd,ntypd,jspd),
     +          ddn(0:lmaxd,ntypd,jspd),kveclo(nlotot),
     +          ulos(nlod,ntypd,jspd),dulos(nlod,ntypd,jspd),
     +          uulon(nlod,ntypd,jspd),dulon(nlod,ntypd,jspd),
     +          evac(2,jspd,nwdd),ello(nlod,ntypd,jspd),bkpt(3),rk(nvd),
     +          el(0:lmaxd,ntypd,jspd,nwdd),k1(nvd),k2(nvd),k3(nvd))

#ifdef CPP_HDF
      CALL openeig(
     >             nbasfcn,neigd,nkpt(1),jspd,lmaxd,nlod,ntypd,
     >             .false.,eonly)
#else
      OPEN (66,file=fname(1),access='direct',form='unformatted',
     +                                     recl=irecl,status='old')
#endif
      IF ( eonly .AND. (irank==0) ) THEN
        OPEN(67,file='eig.soc',form='unformatted',status='replace')
      ENDIF 
c
c---> set up and solve the eigenvalue problem
c
      nk0 = 0
      nk00 = 0
c---> loop over energy windows
      DO 30 nw = 1,nwd
c
c--->    radial k-idp s-o matrix elements calc. and storage
c
         nrec = nk0+1
#ifdef CPP_HDF
         DO jspin = 1, jspins
         CALL read_kptenpar(
     >                      1,jspin,
     <                      bkpt,wk,el(:,:,jspin,nw),ello(:,:,jspin),
     <                      evac(:,jspin,nw))
         ENDDO
#else
         READ (66,rec=nrec) el(:,:,1,nw),evac(:,1,nw),ello(:,:,1)
         IF (jspins.EQ.2) THEN
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
           nrec = nk0+1+nkpt(nw)/isize
           IF (irank.EQ.0) nrec = nk0+1+nkpt(nw)
#else
           nrec = nk0+1+nkpt(nw)
#endif
           READ (66,rec=nrec) el(:,:,jspins,nw),evac(:,jspins,nw),
     +                        ello(:,:,jspins)  
         ENDIF
#endif
         CALL cpu_time(t1)
C  ..
         ALLOCATE( rsopdp(ntypd,lmaxd,2,2),rsopdpd(ntypd,lmaxd,2,2),
     +               rsopp(ntypd,lmaxd,2,2),rsoppd(ntypd,lmaxd,2,2),
     +               rsoplop(ntypd,nlod,2,2),rsoplopd(ntypd,nlod,2,2),
     +               rsopdplo(ntypd,nlod,2,2),rsopplo(ntypd,nlod,2,2),
     +               rsoploplop(ntypd,nlod,nlod,2,2),
     +               soangl(lmaxd,-lmaxd:lmaxd,2,lmaxd,-lmaxd:lmaxd,2) )

         soangl(:,:,:,:,:,:) = cmplx(0.0,0.0)
         CALL spnorb(
     >         ntypd,lmaxd,jmtd,jspd,nwdd,nlhd,nlod,
     >         theta,phi,jspins,ntype,nw,irank,
     >         jri,lmax,dx,rmsh,el,ello,nlo,llo,l_dulo,ulo_der,vr,
     >         soc_opt(ntype+2), 
     <         rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     <         rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <         us,dus,uds,duds,ulos,dulos,uulon,dulon,soangl)
c
         l_all = .false.
         INQUIRE (file='allbut',exist=l_all)
         IF (l_all) THEN
           OPEN (1,file='allbut',form='formatted')
           READ (1,*) n
           WRITE (*,*) 'allbut',n
           CLOSE (1)
           rsopp(1:n-1,:,:,:) = 0.0 ; rsopp(n+1:ntypd,:,:,:) = 0.0 
           rsopdp(1:n-1,:,:,:) = 0.0 ; rsopdp(n+1:ntypd,:,:,:) = 0.0 
           rsoppd(1:n-1,:,:,:) = 0.0 ; rsoppd(n+1:ntypd,:,:,:) = 0.0 
           rsopdpd(1:n-1,:,:,:) = 0.0 ; rsopdpd(n+1:ntypd,:,:,:) = 0.0 
           rsoplop(1:n-1,:,:,:) = 0.0 ; rsoplop(n+1:ntypd,:,:,:) = 0.0 
           rsoplopd(1:n-1,:,:,:) = 0.0 ; rsoplopd(n+1:ntypd,:,:,:) = 0.0 
           rsopdplo(1:n-1,:,:,:) = 0.0 ; rsopdplo(n+1:ntypd,:,:,:) = 0.0 
           rsopplo(1:n-1,:,:,:) = 0.0 ; rsopplo(n+1:ntypd,:,:,:) = 0.0
           rsoploplop(1:n-1,:,:,:,:) = 0.0 
           rsoploplop(n+1:ntypd,:,:,:,:) = 0.0
         ENDIF
         l_all = .false.
         INQUIRE (file='BiSb',exist=l_all)
         IF (l_all) THEN
           OPEN (1,file='BiSb',form='formatted')
           READ (1,*) n
           WRITE (*,*) 'BiSb',n
           CLOSE (1)
           rsopp(:,1,:,:) = 0.0329/2 + 0.0707 * real(n) / 200.0
           rsopp(:,2,:,:) = 0.0077/2 + 0.0265 * real(n) / 200.0
           rsopp(:,3:lmaxd,:,:) = 0.0000 + 0.0005 * real(n) / 200.0
           rsopdp(:,1,:,:) = 0.0125/2 + 0.0178 * real(n) / 200.0
           rsopdp(:,2,:,:) =-0.0009/2 - 0.0075 * real(n) / 200.0
           rsopdp(:,3:lmaxd,:,:) = 0.0000 + 0.0003 * real(n) / 200.0
           rsoppd(:,:,:,:) = rsopdp(:,:,:,:)
           rsopdpd(:,1,:,:) = 0.0046/2 + 0.0043 * real(n) / 200.0
           rsopdpd(:,2,:,:) = 0.0001/2 + 0.0020 * real(n) / 200.0
           rsopdpd(:,3:lmaxd,:,:) = 0.0000 + 0.0002 * real(n) / 200.0
        DO n = 1,ntype
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,1),l=1,3)
         WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,1),l=1,3)
         WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,1),l=1,3)
         WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,1),l=1,3)
        ENDDO
 8000 FORMAT (' spin - orbit parameter HR  ')
 8001 FORMAT (8f8.4)
 9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')

         ENDIF


         IF (irank==0) THEN
           IF (soc_opt(ntype+1) .or. l_all) THEN
             IF (l_all) THEN
               WRITE (6,fmt='(A)') 'Only SOC contribution of certain'
     &          //' atom types included in Hamiltonian.'
             ELSE 
               WRITE (chntype,'(i3)') ntype
               WRITE (6,fmt='(A,2x,'//chntype//'l1)') 'SOC contributi'
     &          //'on of certain atom types included in Hamiltonian:',
     &          (soc_opt(n),n=1,ntype)
             ENDIF
           ELSE
             WRITE(6,fmt='(A,1x,A)') 'SOC contribution of all atom'//
     &        ' types inculded in Hamiltonian.'
           ENDIF
           IF (soc_opt(ntype+2)) THEN
             WRITE(6,fmt='(A)')
     &        'SOC Hamiltonian is constructed by neglecting B_xc.'
           ENDIF
         ENDIF

         ALLOCATE( zso(nbasfcn,2*neigd,jspd),eig_so(2*neigd) )
         zso(:,:,:) = cmplx(0.0,0.0)
         soangl(:,:,:,:,:,:) = conjg(soangl(:,:,:,:,:,:))
         CALL cpu_time(t2)
         cp_spnorb = cp_spnorb + t2 - t1
c
c--->    loop over k-points: each can be a separate task
c
         n_loc = int(nkpt(nw)/isize)
         n_plus = nkpt(nw) - isize*n_loc
         i_plus = -1
         IF (irank.LT.n_plus) i_plus = 0
         n_end = (irank+1)+(n_loc+i_plus)*isize
         nk1 = nk0 
c
c--->  start loop k-pts
c
         DO 20 nk = irank+1,n_end,isize
c
#if ( defined(CPP_MPI) && !defined(CPP_T90) && !defined(CPP_HDF) )
            nk1 = nk1 + 1
            IF (irank.EQ.0) nk1 = nk0 + nk
#else
            nk1 = nk0 + nk
#endif

            bkpt(:) = 0.0 ; rk(:) = 0.0  ! calculate later
            CALL cpu_time(t1)
            CALL alineso(
     >        irank,isize,jspd,nvd,neigd,ntypd,natd,lmaxd,nop,nlod,
     >        nbasfcn,jspins,ntype,soc_opt,neq,lmax,ngopr,invtab,
     >        bmat,omtil,taual,mrot,rmt,
     X        rk,bkpt,
     >        nlo,llo,llod,
     >        rsopp,rsoppd,rsopdp,rsopdpd,nk1,nkpt(nw),l_dulo,
     >        rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     >        us,dus,uds,duds,ulos,dulos,uulon,dulon,soangl,
     >        invsat,invsatnr,ddn,lapw_l,bbmat,l_noco,alph,beta,eonly,
     >        nlotot,odi,ods,
     X        cp_help,cp_ham,cp_diag,
     <        k1,k2,k3,kveclo,nv,wtkpt,ello,nmat,nsz,
     <        eig_so,zso)
            CALL cpu_time(t2)
            cp_aline = cp_aline + t2 - t1
            IF (irank.EQ.0) THEN
              WRITE ( 6,FMT=8010) (bkpt(i),i=1,3),nv,nsz
              WRITE (16,FMT=8010) (bkpt(i),i=1,3),nv,nsz
              WRITE ( 6,FMT=8020) (eig_so(i),i=1,nsz)
              WRITE (16,FMT=8020) (eig_so(i),i=1,nsz)
            ENDIF
 8010 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,' 3-d basis functions',/,
     +       ' the',i4,' SOC eigenvalues are:')
 8020 FORMAT (5x,5f12.6)

            IF (eonly) THEN

              ALLOCATE ( impi(2), rmpi(2*neigd+4) )
              impi(1)= nsz
              impi(2)= nk
              DO j=1,nsz
                rmpi(j)= eig_so(j)
              ENDDO
              rmpi(nsz+1)= wtkpt
              DO j=1,3
                rmpi(nsz+1+j)= bkpt(j)
              ENDDO
#ifdef CPP_MPI
              IF (irank/=0) THEN
                CALL MPI_SEND(impi,2,MPI_INTEGER,0,
     &                        0,MPI_COMM_WORLD,ierr)
                CALL MPI_SEND(rmpi(1),impi(1)+4,TYP_REAL,0,
     &                        0,MPI_COMM_WORLD,ierr)
              ELSE
                DO i = 0,isize-1
                  IF (i/=0) THEN
                    CALL MPI_RECV(impi,2,MPI_INTEGER,i,
     &                            0,MPI_COMM_WORLD,stt,ierr)
                    CALL MPI_RECV(rmpi(1),impi(1)+4,TYP_REAL,i,
     &                            0,MPI_COMM_WORLD,stt,ierr)
                  ENDIF
#endif
                  WRITE(67) impi(2) ! only for test purposes
                  WRITE(67)
     &             rmpi(impi(1)+2), rmpi(impi(1)+3), rmpi(impi(1)+4)
                  WRITE(67) rmpi(impi(1)+1)
                  WRITE(67) impi(1)
                  DO j= 1,impi(1)
                    WRITE(67) rmpi(j)
                  ENDDO
#ifdef CPP_MPI
                ENDDO
              ENDIF
#endif
              DEALLOCATE( impi, rmpi )

            ELSE

              DO jspin = 1,jspins

#if ( defined(CPP_MPI) && !defined(CPP_T90) && !defined(CPP_HDF) )
                IF (irank.NE.0) THEN
                  nrec = nk1
                  IF (jspin.EQ.2) nrec = nk1 + nkpt(nw)/isize
                ELSE
                  nrec = nk0+nk
                  IF (jspin.EQ.2) nrec = nk0+nk+nkpt(nw)
                ENDIF
                nv1(1) = nv ; nv1(2) = nv
               CALL mpi_col_eig(
     >                       irank,isize,jspd,2*neigd,lmaxd,ntypd,nlod,
     >      .false.,l_noco,.false.,jspin,.false.,.false.,film,nkpt(nw),
     X                  nrec,el(:,:,:,nw),evac,ello,bkpt,wtkpt,nsz,nv1,
     X                      nmat,eig_so)
#else
                nrec = nk0+nk
                IF (jspin.EQ.2) nrec = nk0+nk+nkpt(nw)
#endif
#ifdef CPP_HDF
                CALL cpu_time(t1)
                CALL writeeigc(
     >                       nk1,jspin,nsz,nmat,
     >                       eig_so(:nsz),zso(:nmat,:nsz,jspin))
                CALL cpu_time(t2)
                cp_spnorb = t2 - t1
#else
                WRITE (66,rec=nrec) el(:,:,jspin,nw),evac(:,jspin,nw),
     +                  ello(:,:,jspin),bkpt,wtkpt,nsz,nv,nmat,
!-gb     +              eig_so(:),k1,k2,k3,kveclo,conjg(zso(:,:,jspin))
     +                  eig_so(:),k1,k2,k3,kveclo,(zso(:,:,jspin))
#endif
                ! write(*,*) 'process',irank,' wrote ',nrec
              ENDDO

            ENDIF ! (eonly) ELSE

   20     CONTINUE ! DO nk 
          DEALLOCATE (zso,eig_so,rsoploplop,rsopplo,rsopdplo,rsoplopd)
          DEALLOCATE (rsoplop,rsopdp,rsopdpd,rsopp,rsoppd,soangl)
          nk0 = nk0 + jspins*nkpt(nw)
          nk00 = nk00+ nkpt(nw)
   30 CONTINUE

      DEALLOCATE ( vr,us,dus,uds,duds,ulos,dulos,uulon,dulon,ddn )
      DEALLOCATE ( bkpt,rk,k1,k2,k3,evac,ello,el,kveclo )

#ifdef CPP_HDF
      CALL closeeig()
#else
      CLOSE (66)
#endif
      IF ( eonly .AND. (irank==0) ) THEN
        CLOSE(67)
      ENDIF

      IF (irank.EQ.0) THEN
        CALL outtime('SOC:   spnorb :',cp_spnorb)
        CALL outtime('SOC:     -help:',cp_help  )
        CALL outtime('SOC:     -ham :',cp_ham   )
        CALL outtime('SOC:     -diag:',cp_diag  )
        CALL outtime('SOC:   alineso:',cp_aline )
      ENDIF

      RETURN
      END SUBROUTINE eigenso
      END MODULE m_eigenso
