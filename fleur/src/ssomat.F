      SUBROUTINE ssomat(
     > jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,nlod,
     > neigd,nvd,lmaxd,nop,natd,ntypsd,llod,l_dulo,
     > ulo_der,invsat,invsatnr,lapw_l,bbmat,alph,beta,
     > nbasfcn,film,invs,invs2,nq2,nq3,nvac,jspins,
     > nwd,nu,irank,isize,irecl,ntype,soc_opt,nlo,llo,
     > theta,phi,omtil,neq,lmax,mrot,ngopr,nkpt,nlh,
     > jri,ntypsy,bmat,rmt,dx,rmsh,taual,invtab,
     > odi,ods,qss,zelec,tkb,sso_opt,TYP_REAL )
c
#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : openeig, closeeig, read_kpt, read_enpar, 
     &                        read_eig, read_vec, read_basis
#endif
      USE m_spnorb 
      USE m_loddop
      USE m_abcof 
      USE m_od_types, ONLY : od_inp, od_sym
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd
      INTEGER, INTENT (IN) :: neigd,nvd,lmaxd,nop,natd,llod
      INTEGER, INTENT (IN) :: ntype,nu,irecl,nq2,nq3,nvac,jspins
      INTEGER, INTENT (IN) :: irank,isize,ntypsd,nlod,nwd 
      INTEGER, INTENT (IN) :: nbasfcn
      REAL,    INTENT (IN) :: theta,phi,omtil
      LOGICAL, INTENT (IN) :: film,invs,invs2
      INTEGER, INTENT (IN) :: TYP_REAL 
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),nkpt
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),jri(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd),lapw_l(ntype)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd)
      REAL,    INTENT (IN) :: bmat(3,3),dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: rmt(ntypd),taual(3,natd),bbmat(3,3)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd)
      LOGICAL, INTENT (IN) :: soc_opt(ntypd+2),sso_opt(2) 
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      REAL,    INTENT (IN) :: qss(3)
      REAL,    INTENT (IN) :: zelec, tkb  
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
C     ..
C     .. Locals ..
      LOGICAL :: eachtype,twosteps, merging, matform, step1, step2   
      INTEGER :: matelversion 

      INTEGER :: lmd 
      INTEGER :: n_loc,n_plus,i_plus,n_end,nnode,  
     &           ne,nmat,nlotot,nv(2)
      INTEGER :: i1,i2,js,jsloc,nk,n,l,m,lm,band,cyc,cycnode,kp
      INTEGER :: cyclen,ncyc,ncycnode 
      INTEGER :: absame,abnext    
      INTEGER :: nat,na,natall 
      REAL    :: r1, r2, pi, bkpt(3), wtkpt  
      REAL    :: eps 
      COMPLEX :: c1 
      CHARACTER(len=8) :: dop,iop,name(10)
      CHARACTER(len=8) :: abfile 
      CHARACTER(len=3) :: chntype

      REAL,    ALLOCATABLE :: bkp(:,:)
      INTEGER, ALLOCATABLE :: knext(:), cycstart(:)
      COMPLEX, ALLOCATABLE :: matel(:,:,:) 


      REAL,    ALLOCATABLE :: eig(:,:,:), eignode(:,:,:), eigk(:) 
      INTEGER, ALLOCATABLE :: knr(:,:), knrnode(:,:)   
      REAL,    ALLOCATABLE :: wkpt(:), wkptnode(:) 
      COMPLEX, ALLOCATABLE :: soangl(:,:,:,:,:,:)
      REAL,    ALLOCATABLE :: rsopp(:,:,:,:),rsoppd(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopdp(:,:,:,:),rsopdpd(:,:,:,:)
      COMPLEX, ALLOCATABLE :: acof(:,:,:,:,:,:), bcof(:,:,:,:,:,:)

      REAL,    ALLOCATABLE :: vr(:,:,:,:)
      REAL,    ALLOCATABLE :: us(:,:,:),dus(:,:,:)
      REAL,    ALLOCATABLE :: uds(:,:,:),duds(:,:,:)
      REAL,    ALLOCATABLE :: ddn(:,:,:) 
      REAL,    ALLOCATABLE :: evac(:,:),ello(:,:,:)
      REAL,    ALLOCATABLE :: el(:,:,:),elspnorb(:,:,:,:)
      INTEGER, ALLOCATABLE :: k1(:,:),k2(:,:),k3(:,:)

      REAL,    ALLOCATABLE :: vz(:,:,:)
      COMPLEX, ALLOCATABLE :: vzxy(:,:,:,:),vpw(:,:)
      COMPLEX, ALLOCATABLE :: ccof(:,:,:,:) 
      INTEGER, ALLOCATABLE :: kveclo(:)
      REAL,    ALLOCATABLE :: ulos(:,:,:),dulos(:,:,:) 
      REAL,    ALLOCATABLE :: uulon(:,:,:),dulon(:,:,:)
      REAL,    ALLOCATABLE :: rsoplop(:,:,:,:)
      REAL,    ALLOCATABLE :: rsoplopd(:,:,:,:),rsopdplo(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopplo(:,:,:,:),rsoploplop(:,:,:,:,:)

#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
! this subroutine is never called if CPP_INVERSION,
! but the compiler expects this statement:  
      REAL, ALLOCATABLE :: z(:,:)
#else
      COMPLEX, ALLOCATABLE :: z(:,:)
#endif 

#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER ierr(3), stt(MPI_STATUS_SIZE)
#endif

C     ..................................................................

      merging= .false. 
      ! merging= .true. => all matrix-element-files will be merged into
      !       one big file (in case of huge files that takes much time)  

      matform= .false. 
      ! matform=.true. => file of matrix-elements formatted  


      matelversion= 1 

      eps= 1.e-15 

      pi= 4.*ATAN(1.) 

      lmd= lmaxd*(lmaxd+2)
      nlotot= 0
      DO n= 1,ntype
        DO l= 1,nlo(n)
          nlotot= nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        ENDDO
      ENDDO

      IF (nlotot /= 0) THEN
        STOP 'Stop in ssomat:  local orbitals not implemented'
      ENDIF 
      IF (nwd /= 1) THEN
        STOP 'Stop in ssomat:  more than one window never tested'   
      ENDIF 

      twosteps= sso_opt(1) 
      eachtype= sso_opt(2)

      natall= 0 
      DO n= 1,ntype 
        natall= natall+neq(n)
      ENDDO 

      IF (twosteps) THEN
        INQUIRE(file='abcof000',exist=step2)
        step1= (.not. step2)
      ELSE
        step1= .true.
        step2= .true.
      ENDIF 


#ifdef CPP_HDF
      CALL openeig(
     > nbasfcn,neigd,nkpt,jspd,lmaxd,nlod,ntypd,
     > .false.,.true. )
#else
      OPEN(66,file='eig',
     & access='direct',form='unformatted',recl=irecl,
     & status='old',action='read')
#endif
      OPEN (nu,file='pottot',form='unformatted',
     & status='old',action='read')


      ALLOCATE( bkp(3,nkpt), knext(nkpt), cycstart(nkpt) )


      IF (irank==0) THEN
        WRITE(*,fmt='(A)')
     &   'reading position of k-points in eig-file ...'
#ifdef CPP_HDF 
        DO nk= 1,nkpt
          CALL read_kpt(nk,1, bkpt,wtkpt)
          bkp(1,nk)= bkpt(1) ; bkp(2,nk)=bkpt(2) ; bkp(3,nk)= bkpt(3)
        ENDDO
#else
        ALLOCATE( evac(2,jspd),ello(nlod,ntypd,jspd),
     &            el(0:lmaxd,ntypd,jspd) )
        DO nk= 1,nkpt
          READ(66,rec=nk) el,evac,ello, bkpt
          bkp(1,nk)= bkpt(1) ; bkp(2,nk)=bkpt(2) ; bkp(3,nk)= bkpt(3)
        ENDDO
        DEALLOCATE( evac,ello,el )
#endif 
        WRITE(*,fmt='(A)') ' ... done'
      ENDIF
#ifdef CPP_MPI
      CALL MPI_BCAST(bkp,3*nkpt,TYP_REAL,0,MPI_COMM_WORLD,ierr)
#endif

      CALL nextk(nkpt,bkp,qss, knext,ncyc,cyclen,cycstart)

      IF (irank<ncyc) THEN 


      ncycnode= ncyc/isize
      IF ( MOD(ncyc,isize) > irank ) THEN
        ncycnode= ncycnode + 1
      ENDIF

      ALLOCATE( eigk(neigd) ) 

      ! needed directly for calculating matrix elements  

      ALLOCATE( soangl(lmaxd,-lmaxd:lmaxd,2,lmaxd,-lmaxd:lmaxd,2) )
      ALLOCATE( rsopdp(ntypd,lmaxd,2,2),rsopdpd(ntypd,lmaxd,2,2),
     +          rsopp(ntypd,lmaxd,2,2),rsoppd(ntypd,lmaxd,2,2) )
      ALLOCATE( acof(neigd,0:lmd,natd,2,2,2),
     +          bcof(neigd,0:lmd,natd,2,2,2) )
                !   (neigd, 0:lmd, natd , local-spin, global-spin, kp) 

      ! needed to calculate a,bcofs or k-independent integrals 

      ALLOCATE( vr(jmtd,0:nlhd,ntypd,jspd) ) 
      ALLOCATE(  us(0:lmaxd,ntypd,jspd), dus(0:lmaxd,ntypd,jspd),
     +          uds(0:lmaxd,ntypd,jspd),duds(0:lmaxd,ntypd,jspd),
     +          ddn(0:lmaxd,ntypd,jspd) )
      ALLOCATE( z(nbasfcn,neigd) ) 
      ALLOCATE( evac(2,jspd),ello(nlod,ntypd,jspd),
     +          el(0:lmaxd,ntypd,jspd),
     +          elspnorb(0:lmaxd,ntypd,jspd,1),
     +          k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd) )

      ! dummies (mainly for lo)  

      ALLOCATE( vz(nmzd,2,jspd),
     +          vzxy(nmzxyd,n2d-1,2,jspd),vpw(nq3,jspd) )
      ALLOCATE( ccof(-llod:llod,neigd,nlod,natd) )
      ALLOCATE( kveclo(nlotot) )
      ALLOCATE( ulos(nlod,ntypd,jspd),dulos(nlod,ntypd,jspd),
     +          uulon(nlod,ntypd,jspd),dulon(nlod,ntypd,jspd) )
      ALLOCATE( rsoplop(ntypd,nlod,2,2),rsoplopd(ntypd,nlod,2,2),
     +          rsopdplo(ntypd,nlod,2,2),rsopplo(ntypd,nlod,2,2),
     +          rsoploplop(ntypd,nlod,nlod,2,2) )

      CALL loddop(
     > jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     > jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     > nlh,jri,ntypsd,ntypsy,nu,natd,neq,
     < iop,dop,i1,vr,vpw,vz,vzxy,name)


#ifdef CPP_HDF
      CALL read_enpar(1,1, el(:,:,1),ello(:,:,1),evac(:,1)) 
      CALL read_enpar(1,2, el(:,:,2),ello(:,:,2),evac(:,2)) 
#else
      READ (66,rec=1) el,evac,ello
#endif 

      ! elspnorb probably not necessary !
      DO l= 0,lmaxd
        DO n= 1,ntypd
          DO jsloc= 1,jspd
            elspnorb(l,n,jsloc,1)= el(l,n,jsloc)
          ENDDO
        ENDDO
      ENDDO 
      ! use same definition of rotation angles as in noco-routines 
      r1= -theta
      r2= phi+pi
      CALL spnorb(
     > ntypd,lmaxd,jmtd,jspd,1,nlhd,nlod,
     > r1,r2,jspins,ntype,1,irank,
     > jri,lmax,dx,rmsh,elspnorb,ello,nlo,llo,l_dulo,
     > ulo_der,vr,.true., 
     < rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     < rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     < us,dus,uds,duds,ulos,dulos,uulon,dulon,soangl)
      IF ( (irank==0) .and. step2 ) THEN 
        IF (soc_opt(ntype+1)) THEN
          WRITE (chntype,'(i3)') ntype
          WRITE (6,fmt='(A,2x,'//chntype//'l1)') 
     &     'SOC contribution of certain atom types included in'// 
     &     ' Hamiltonian:', (soc_opt(n),n=1,ntype)
        ELSE
          WRITE(6,fmt='(A)') 
     &     'SOC contribution of all atom types inculded in Hamiltonian.'
        ENDIF
        IF (.not. soc_opt(ntype+2)) THEN
          WRITE(6,fmt='(A)') 'For ss+soc, spav is always set to true.'
        ENDIF 
      ENDIF 

      IF (step1) THEN 

      abfile= 'abcofXXX'
      abfile(6:6)= ACHAR(IACHAR('0')+irank/100 )
      abfile(7:7)= ACHAR(IACHAR('0')+MOD(irank,100)/10 )
      abfile(8:8)= ACHAR(IACHAR('0')+MOD(irank,10) )
      OPEN(301,file=abfile,form='unformatted',status='replace')
      WRITE(301) isize 

      DO cycnode= 1,ncycnode 

        cyc= (cycnode-1)*isize+irank+1 

        nk= cycstart(cyc)
        DO kp= 1,cyclen+1

#ifdef CPP_HDF
          CALL read_eig(nk,1, ne) 
          CALL read_kpt(nk,1, bkpt,wtkpt)
#if ( defined(CPP_INVERSION) )
          STOP 'ssomat: Compile without inversion!'
#else
          CALL read_vec(nk,1,ne,1,ne,nmat,z)
#endif
          CALL read_basis(nk,1,nlotot, 
     <     nmat,nv(1),k1(:,1),k2(:,1),k3(:,1),kveclo)
          CALL read_basis(nk,2,nlotot, 
     <     nmat,nv(2),k1(:,2),k2(:,2),k3(:,2),kveclo)
#else
          READ (66,rec=nk) el,evac,ello,
     &                     bkpt,wtkpt,ne,nv,nmat,
     &                     eigk,k1,k2,k3,z
#endif 
          IF (ne/=neigd) THEN
            STOP 'Stop in ssomat:  ne /= neigd'
          ENDIF

          DO jsloc= 1,2 
            CALL abcof(
     >       lmaxd,ntypd,neigd,neigd,natd,nop,nvd,jspd,
     >       lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >       ntype,mrot,ngopr,taual,neq,lmax,rmt,omtil,bmat,
     >       bbmat,bkpt,k1,k2,k3,nv,nmat,neigd,z,
     >       us(0,1,jsloc),dus(0,1,jsloc),uds(0,1,jsloc),
     >       duds(0,1,jsloc),ddn(0,1,jsloc),invsat,invsatnr,
     >       ulos(1,1,jsloc),uulon(1,1,jsloc),dulon(1,1,jsloc),
     >       dulos(1,1,jsloc),llo,nlo,l_dulo,lapw_l,
     >       .true.,.true.,jsloc,alph,beta,qss,kveclo,odi,ods,
     <       acof(1,0,1,jsloc,1,2),bcof(1,0,1,jsloc,1,2),ccof)
          ENDDO

          ! rotate abcof(...,2) in global spin coordinate frame, 
          ! write result to abcof(...,1)
          nat= 0 
          DO n= 1,ntype  
            DO na= 1,neq(n)
              nat= nat+1 
              r1= alph(n) 
              r2= beta(n) 
              DO lm= 0,lmd
                DO band= 1,neigd
                  acof(band,lm,nat,1,1,1)= 
     &                CMPLX( COS(r1/2.) , -SIN(r1/2.) )
     &               *CMPLX( COS(r2/2.) , 0. )  
     &               *acof(band,lm,nat,1,1,2)
                  acof(band,lm,nat,2,1,1)= 
     &                CMPLX( COS(r1/2.) , -SIN(r1/2.) )
     &               *CMPLX( -SIN(r2/2.) , 0. )  
     &               *acof(band,lm,nat,2,1,2)
                  acof(band,lm,nat,1,2,1)= 
     &                CMPLX( COS(r1/2.) , +SIN(r1/2.) )
     &               *CMPLX( +SIN(r2/2.) , 0. )  
     &               *acof(band,lm,nat,1,1,2)
                  acof(band,lm,nat,2,2,1)= 
     &                CMPLX( COS(r1/2.) , +SIN(r1/2.) )
     &               *CMPLX( COS(r2/2.) , 0. )  
     &               *acof(band,lm,nat,2,1,2)
                  bcof(band,lm,nat,1,1,1)= 
     &                CMPLX( COS(r1/2.) , -SIN(r1/2.) )
     &               *CMPLX( COS(r2/2.) , 0. )  
     &               *bcof(band,lm,nat,1,1,2)
                  bcof(band,lm,nat,2,1,1)= 
     &                CMPLX( COS(r1/2.) , -SIN(r1/2.) )
     &               *CMPLX( -SIN(r2/2.) , 0. )  
     &               *bcof(band,lm,nat,2,1,2)
                  bcof(band,lm,nat,1,2,1)= 
     &                CMPLX( COS(r1/2.) , +SIN(r1/2.) )
     &               *CMPLX( +SIN(r2/2.) , 0. )  
     &               *bcof(band,lm,nat,1,1,2)
                  bcof(band,lm,nat,2,2,1)= 
     &                CMPLX( COS(r1/2.) , +SIN(r1/2.) )
     &               *CMPLX( COS(r2/2.) , 0. )  
     &               *bcof(band,lm,nat,2,1,2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          WRITE(301) nk
          WRITE(301) bkpt(1), bkpt(2), bkpt(3), wtkpt
          ! acof(neigd,0:lmd,natall,2,2,2)
          DO js= 1,2
            DO jsloc= 1,2
              DO n= 1,natall 
                DO lm= 0,lmd
                  DO band= 1,neigd  
                    WRITE(301) acof(band,lm,n,jsloc,js,1),
     &                         bcof(band,lm,n,jsloc,js,1) 
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO 

          nk= knext(nk)
        ENDDO ! kp= 1,cyclen+1 

      ENDDO ! cycnode 

      CLOSE(301)

      ENDIF ! (step1)

      DEALLOCATE( 
     & vr, us,dus,uds,duds,ddn, z, 
     & evac,ello,el,elspnorb,k1,k2,k3 )
      DEALLOCATE( 
     & vz,vzxy,vpw, ccof, kveclo, ulos,dulos,uulon,dulon, 
     & rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop )

      IF (step2) THEN 

      ALLOCATE( eignode(neigd,cyclen,ncycnode), 
     &          knrnode(cyclen,ncycnode), wkptnode(ncycnode) )

      abfile= 'abcofXXX'
      abfile(6:6)= ACHAR(IACHAR('0')+irank/100 )
      abfile(7:7)= ACHAR(IACHAR('0')+MOD(irank,100)/10 )
      abfile(8:8)= ACHAR(IACHAR('0')+MOD(irank,10) )
      OPEN(301,file=abfile,form='unformatted') 
      READ(301) i1
      IF (i1/=isize) THEN
        CLOSE(301)
        WRITE (*,fmt='(A,x,A,i3,A,i3)')
     &   abfile,'was written when isize =',i1,', now isize =',isize
        STOP 'Stop in ssomat'
      ENDIF

      abfile= 'matelXXX'
      abfile(6:6)= ACHAR(IACHAR('0')+irank/100 )
      abfile(7:7)= ACHAR(IACHAR('0')+MOD(irank,100)/10 )
      abfile(8:8)= ACHAR(IACHAR('0')+MOD(irank,10) )
      OPEN(302,file=abfile,
     & form='unformatted',status='replace') 
      ALLOCATE( matel(neigd,neigd,0:ntype) )
#ifndef CPP_HDF
      ALLOCATE( evac(2,jspd),ello(nlod,ntypd,jspd),
     &          el(0:lmaxd,ntypd,jspd) )
#endif 

      DO cycnode= 1,ncycnode 

        cyc= (cycnode-1)*isize+irank+1

        WRITE(*,fmt='(A,i4,x,A,i4)') 
     &   'matrix elements of cycle',cyc,'/',ncyc 
        READ(301) nk
        READ(301) bkpt(1), bkpt(2), bkpt(3), wtkpt
        DO js= 1,2
          DO jsloc= 1,2
            DO n= 1,natall 
              DO lm= 0,lmd
                DO band= 1,neigd
                  READ(301) acof(band,lm,n,jsloc,js,1),
     &                      bcof(band,lm,n,jsloc,js,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO 

        wkptnode(cycnode)= wtkpt  

        DO kp= 1,cyclen 
          absame= MOD(kp+1,2)+1
          abnext= MOD(kp,2)+1

#ifdef CPP_HDF
          CALL read_eig(nk,1, ne,eigk) 
#else
          READ (66,rec=nk) el,evac,ello,
     &                     bkpt,r1,ne,nv,nmat,
     &                     eigk 
#endif
          DO band= 1,neigd
            eignode(band,kp,cycnode)= eigk(band)
          ENDDO 

          IF (ne/=neigd) THEN
            STOP 'Stop in ssomat:  ne /= neigd'
          ENDIF
          knrnode(kp,cycnode)= nk 

          WRITE(302) nk, bkpt(1), bkpt(2), bkpt(3)  

          ! matrix elements within k
          CALL ssomatel(
     >     neigd,lmaxd,lmd,ntypd,natd, 
     >     lmax,ntype,neq,soc_opt, 
     >     soangl,rsopp,rsoppd,rsopdp,rsopdpd, 
     >     .true.,
     >     acof(1,0,1,1,1,absame),bcof(1,0,1,1,1,absame),
     >     acof(1,0,1,1,1,absame),bcof(1,0,1,1,1,absame),
     <     matel )
          DO n= 0,ntype 
            IF ( eachtype .or. (n==0) ) THEN
              DO band= 1,neigd
                DO i1= 1,band  
                  WRITE(302) matel(i1,band,n)    
                ENDDO 
              ENDDO
            ENDIF
          ENDDO

          READ(301) nk
          READ(301) bkpt(1), bkpt(2), bkpt(3), r1 
          IF ( ABS(r1-wtkpt) > eps ) THEN
            STOP 'Stop in ssomat:  weights of k-points differ'
          ENDIF
          DO js= 1,2
            DO jsloc= 1,2
              DO n= 1,natall 
                DO lm= 0,lmd
                  DO band= 1,neigd
                    READ(301) acof(band,lm,n,jsloc,js,abnext),
     &                        bcof(band,lm,n,jsloc,js,abnext)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          ! matrix elements between k and k+q 
          CALL ssomatel(
     >     neigd,lmaxd,lmd,ntypd,natd, 
     >     lmax,ntype,neq,soc_opt, 
     >     soangl,rsopp,rsoppd,rsopdp,rsopdpd, 
     >     .false., 
     >     acof(1,0,1,1,1,absame),bcof(1,0,1,1,1,absame),
     >     acof(1,0,1,1,1,abnext),bcof(1,0,1,1,1,abnext),
     <     matel )
          DO n= 0,ntype
            IF ( eachtype .or. (n==0) ) THEN
              DO band= 1,neigd
                DO i1= 1,neigd 
                  WRITE(302) matel(i1,band,n)    
                ENDDO 
              ENDDO
            ENDIF
          ENDDO

        ENDDO ! kp= 1,cyclen

      ENDDO ! cycnode= 1,ncycnode  

      IF (step1) THEN 
        CLOSE(301,status='delete') 
      ELSE 
        CLOSE(301) 
      ENDIF 
      CLOSE(302) 

      DEALLOCATE( matel )
#ifndef CPP_HDF
      DEALLOCATE( evac,ello,el )
#endif 

      IF (irank==0) THEN 
        ALLOCATE( eig(neigd,cyclen,ncyc), 
     &            knr(cyclen,ncyc), wkpt(ncyc) )
      ENDIF 

#ifdef CPP_MPI
      IF (irank/=0) THEN
        DO cycnode= 1,ncycnode
          CALL MPI_SEND(eignode(1,1,cycnode),
     &     neigd*cyclen,TYP_REAL,0,0,MPI_COMM_WORLD,ierr)  
          CALL MPI_SEND(knrnode(1,cycnode),
     &     cyclen,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)  
          CALL MPI_SEND(wkptnode(cycnode),
     &     1,TYP_REAL,0,0,MPI_COMM_WORLD,ierr)  
        ENDDO
      ENDIF
#endif 
      IF (irank==0) THEN
        DO cyc= 1,ncyc
          nnode= MOD(cyc-1,isize)
          IF (nnode==0) THEN
            cycnode= (cyc-1)/isize+1
            wkpt(cyc)= wkptnode(cycnode)
            DO kp= 1,cyclen
              DO band= 1,neigd
                eig(band,kp,cyc)= eignode(band,kp,cycnode)
              ENDDO 
              knr(kp,cyc)= knrnode(kp,cycnode) 
            ENDDO 
#ifdef CPP_MPI
          ELSE
            CALL MPI_RECV(eig(1,1,cyc),
     &       neigd*cyclen,TYP_REAL,nnode,0,MPI_COMM_WORLD,stt,ierr)
            CALL MPI_RECV(knr(1,cyc),
     &       cyclen,MPI_INTEGER,nnode,0,MPI_COMM_WORLD,stt,ierr)
            CALL MPI_RECV(wkpt(cyc),
     &       1,TYP_REAL,nnode,0,MPI_COMM_WORLD,stt,ierr)
#endif
          ENDIF 
        ENDDO 
      ENDIF ! (irank==0) 

      DEALLOCATE( eignode, knrnode, wkptnode ) 

      ENDIF ! (step2)

      DEALLOCATE( eigk, soangl, rsopdp,rsopdpd,rsopp,rsoppd, acof,bcof )


      ENDIF ! (irank<ncyc)

      DEALLOCATE( bkp, knext, cycstart ) 

#ifdef CPP_HDF
      CALL closeeig()
#else 
      CLOSE(66)
#endif 
      CLOSE(nu)


#ifdef CPP_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      IF ( (irank==0) .and. step2 ) THEN

        IF (matform) THEN
          OPEN(300,file='matelsso_form',
     &     form='formatted',status='replace')
          WRITE(300,fmt='(A)') 'ssomat'
          WRITE(300,fmt='(A,i4)') 'version =',matelversion 
          WRITE(300,fmt='(4i8,x,l1)')
     &     ncyc, cyclen, neigd, ntype, eachtype
          WRITE(300,fmt='(3(e24.18,x))') qss(1), qss(2), qss(3)
          WRITE(300,fmt='(2(e24.18,x))') theta/pi, phi/pi
          WRITE(300,fmt='(2(e24.18,x))') zelec, tkb
          WRITE(300,fmt='(l1,i4)') merging, isize  
          DO cyc= 1,ncyc  
            WRITE(300,fmt='(i8,x,e24.18)') cyc, wkpt(cyc) 
            DO kp= 1,cyclen
              WRITE(300,fmt='(A,i8)') 'knr =',knr(kp,cyc)  
              DO band= 1,neigd
                WRITE(300,fmt='(i8,x,e24.18)') band, eig(band,kp,cyc)
              ENDDO
            ENDDO
          ENDDO  
          IF (merging) THEN
            WRITE(300,fmt='(A)') 'matrix-elements:'
          ENDIF  
        ELSE
          OPEN(300,file='matelsso',
     &     form='unformatted',status='replace')
          WRITE(300) 'ssomat'
          WRITE(300) matelversion  
          WRITE(300)
     &     ncyc, cyclen, neigd, ntype, eachtype
          WRITE(300) qss(1), qss(2), qss(3)
          WRITE(300) theta/pi, phi/pi
          WRITE(300) zelec, tkb
          WRITE(300) merging, isize  
          DO cyc= 1,ncyc  
            WRITE(300) wkpt(cyc)  
            DO kp= 1,cyclen
              WRITE(300) knr(kp,cyc)
              DO band= 1,neigd
                WRITE(300) eig(band,kp,cyc)
              ENDDO
            ENDDO
          ENDDO  
        ENDIF

        DEALLOCATE( eig, knr, wkpt ) 

        IF (merging) THEN

          WRITE (*,fmt='(A)') 'merging files ...'

          DO nnode= 0, isize-1
            IF (nnode<ncyc) THEN 
              abfile= 'matelXXX'
              abfile(6:6)= ACHAR(IACHAR('0')+nnode/100 )
              abfile(7:7)= ACHAR(IACHAR('0')+MOD(nnode,100)/10 )
              abfile(8:8)= ACHAR(IACHAR('0')+MOD(nnode,10) )
              OPEN(301+nnode,file=abfile,
     &         form='unformatted',action='read')
            ENDIF 
          ENDDO 

          DO cyc= 1,ncyc
            nnode= MOD(cyc-1,isize) 

            DO kp= 1,cyclen

              READ(301+nnode) nk, bkpt(1), bkpt(2), bkpt(3)
              IF (matform) THEN 
                WRITE(300,fmt='(A,i8,3(x,e24.18))')
     &           'kpoint', nk, bkpt(1), bkpt(2), bkpt(3)
              ELSE
                WRITE(300)
     &           nk, bkpt(1), bkpt(2), bkpt(3)
              ENDIF 

              DO n= 0,ntype  
                IF ( eachtype .or. (n==0) ) THEN 
                  DO band= 1,neigd
                    DO i1= 1,band
                      READ(301+nnode) c1
                      IF (matform) THEN 
                        WRITE(300,fmt='(2(e24.18,x))') 
     &                   REAL(c1),AIMAG(c1)
                      ELSE 
                        WRITE(300) c1 
                      ENDIF
                    ENDDO 
                  ENDDO
                ENDIF
              ENDDO
              DO n= 0,ntype  
                IF ( eachtype .or. (n==0) ) THEN 
                  DO band= 1,neigd
                    DO i1= 1,neigd
                      READ(301+nnode) c1
                      IF (matform) THEN 
                        WRITE(300,fmt='(2(e24.18,x))') 
     &                   REAL(c1),AIMAG(c1)
                      ELSE 
                        WRITE(300) c1 
                      ENDIF
                    ENDDO 
                  ENDDO
                ENDIF
              ENDDO

            ENDDO ! kp= 1,cyclen 

          ENDDO ! cyc= 1,ncyc  

          DO nnode= 0, isize-1
            CLOSE(301+nnode,status='delete')  
          ENDDO 

          WRITE(*,fmt='(A)') ' ... done'

        ENDIF ! (merging) 

        CLOSE(300)

      ENDIF ! ( (irank==0) .and. step2 )

#ifdef CPP_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      IF (irank==0) THEN  
        IF (step2) THEN 
          IF (merging) THEN  
            WRITE(*,fmt='(A)') 'matrix elements written to "matelsso"'
          ELSE 
            WRITE(*,fmt='(A)') '"matel*" written'
          ENDIF 
        ELSE
          WRITE(*,fmt='(A)') 'abcof-files written'
        ENDIF 
      ENDIF 

      END SUBROUTINE ssomat 
 
! ==================================================================== ! 

      SUBROUTINE nextk(nkpt,bk,qss, knext,ncyc,cyclen,cycstart)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: nkpt
      REAL,    INTENT(IN)  :: bk(3,nkpt)
      REAL,    INTENT(IN)  :: qss(3)

      INTEGER, INTENT(OUT) :: knext(nkpt), ncyc, cyclen, cycstart(nkpt) 

      REAL    :: eps
      INTEGER :: kin, kpq, kp, xy
      REAL    :: kipq(3), kt(3)

      INTEGER :: cyc 
      LOGICAL :: kdone(nkpt)


      eps= 1.e-10

      DO kin= 1,nkpt

        IF ( ( ABS(qss(1)) < eps ) .and.  
     &       ( ABS(qss(2)) < eps ) .and.  
     &       ( ABS(qss(3)) < eps )          ) THEN
   
          kpq= kin

        ELSE

          DO xy= 1,3
            kipq(xy)= bk(xy,kin)+qss(xy)
            DO WHILE ( kipq(xy)-eps/2. > +1./2. )
              kipq(xy)= kipq(xy)-1.
            ENDDO
            DO WHILE ( kipq(xy)-eps < -1./2. )
              kipq(xy)= kipq(xy)+1.
            ENDDO
          ENDDO

          kpq= -1

          DO kp= 1,nkpt
 
            DO xy= 1,3
              kt(xy)= bk(xy,kp)
              DO WHILE ( kt(xy)-eps/2. > +1./2. )
                kt(xy)= kt(xy)-1.
              ENDDO
              DO WHILE ( kt(xy)-eps < -1./2. )
                kt(xy)= kt(xy)+1.
              ENDDO
            ENDDO

            IF ( ( ABS(kipq(1)-kt(1)) < eps ) .and. 
     &           ( ABS(kipq(2)-kt(2)) < eps ) .and.  
     &           ( ABS(kipq(3)-kt(3)) < eps )          ) THEN
              IF (kpq<0) THEN
                kpq= kp
              ELSE
                STOP 'Stop in nextk:  identical k-points detected'
              ENDIF
            ENDIF

          ENDDO

        ENDIF 

        IF (kpq<0) THEN
          STOP 'Stop in nextk: k-point grid does not match with qss'
        ENDIF

        knext(kin)= kpq

      ENDDO  ! kin


      DO kp= 1,nkpt
        kdone(kp)= .false.
        cycstart(kp)= 0
      ENDDO 

      cyclen= 1 
      cycstart(1)= 1
      kdone(1)= .true. 
      kp= knext(1)  
      DO WHILE (kp/=1) 
        IF (kdone(kp)) THEN
          STOP 'Stop in nextk: error 1' 
        ENDIF
        kdone(kp)= .true. 
        IF (cyclen == nkpt) THEN
          STOP 'Stop in nextk: error 2' 
        ENDIF 
        cyclen= cyclen+1
        kp= knext(kp)
      ENDDO 
      IF ( MOD(nkpt,cyclen)/=0 ) THEN 
        STOP 'Stop in nextk: error 3' 
      ENDIF
      ncyc= nkpt/cyclen 

      IF (ncyc>1) THEN 
        DO cyc= 2,ncyc
          kp= cycstart(cyc-1)+1
          IF (kp>nkpt) THEN 
            STOP 'Stop in nextk: error 4'
          ENDIF
          DO WHILE (kdone(kp))
            kp= kp+1
            IF (kp>nkpt) THEN 
              STOP 'Stop in nextk: error 5'
            ENDIF
          ENDDO
          cycstart(cyc)= kp
          kdone(kp)= .true.
          IF (cyclen>1) THEN
            DO kpq= 2,cyclen
              kp= knext(kp)
              IF (kdone(kp)) THEN
                STOP 'Stop in nextk: error 6'
              ENDIF
              kdone(kp)= .true.
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      ! check, whether all k-points classified  
      DO kp= 1,nkpt
        IF ( .not. kdone(kp) ) THEN 
          STOP 'Stop in nextk: error 7'
        ENDIF
      ENDDO

      END SUBROUTINE nextk

! ==================================================================== ! 

      SUBROUTINE ssomatel(
     > neigd,lmaxd,lmd,ntypd,natd, 
     > lmax,ntype,neq,soc_opt,
     > soangl,rsopp,rsoppd,rsopdp,rsopdpd, 
     > diag, 
     > acof1,bcof1,acof2,bcof2, 
     < matel )

      IMPLICIT NONE

      LOGICAL, INTENT(IN)  :: diag 
      INTEGER, INTENT(IN)  :: neigd,lmaxd,lmd,ntypd,natd 
      INTEGER, INTENT(IN)  :: lmax(ntypd),ntype,neq(ntypd) 
      LOGICAL, INTENT(IN)  :: soc_opt(ntypd+2) 
      REAL,    INTENT(IN)  :: 
     &          rsopp(ntypd,lmaxd,2,2), rsoppd(ntypd,lmaxd,2,2),
     &          rsopdp(ntypd,lmaxd,2,2), rsopdpd(ntypd,lmaxd,2,2)  
      COMPLEX, INTENT(IN)  :: 
     &          soangl(lmaxd,-lmaxd:lmaxd,2,lmaxd,-lmaxd:lmaxd,2),  
     &          acof1(neigd,0:lmd,natd,2,2), 
     &          bcof1(neigd,0:lmd,natd,2,2),
     &          acof2(neigd,0:lmd,natd,2,2), 
     &          bcof2(neigd,0:lmd,natd,2,2)
                !    (neigd, 0:lmd, natd , local-spin, global-spin)

      COMPLEX, INTENT(OUT) :: matel(neigd,neigd,0:ntype)

      INTEGER :: band1,band2, n,nat,na, l,m1,m2,lm1,lm2,
     &           jsloc1,jsloc2, js1,js2,jsnumber     
      COMPLEX :: sa(2),sb(2), ra(2),rb(2) 

      ! with the following nesting of loops the calculation of the 
      ! matrix-elements is of order
      ! natall*lmd*neigd*(lmd+neigd) ; note that  lmd+neigd << lmd*neigd

      DO n= 0,ntype 
        DO band2= 1,neigd
          DO band1= 1,neigd
            matel(band1,band2,n)= CMPLX(0.,0.) 
          ENDDO
        ENDDO
      ENDDO 

      ! within one k-point loop over global spin 
      IF (diag) THEN 
        jsnumber= 2
      ELSE 
        jsnumber= 1
      ENDIF
      DO js2= 1,jsnumber 
        IF (diag) THEN
          js1= js2
        ELSE
          js1= 2
        ENDIF

        ! loop over MT 
        na= 0 
        DO n= 1,ntype  
          DO nat= 1,neq(n) 
            na= na+1 
            ! loop over l
            DO l= 1,lmax(n)

              ! loop over m1 
              DO m1= -l,l
                lm1= l*(l+1) + m1

                ! loop over eigenstates
                DO band2= 1,neigd

                  DO jsloc2= 1,2
                    sa(jsloc2)= CMPLX(0.,0.)
                    sb(jsloc2)= CMPLX(0.,0.)
                    DO m2= -l,l
                      lm2= l*(l+1) + m2

                      sa(jsloc2)= sa(jsloc2) + 
     &                   CONJG(acof2(band2,lm2,na,jsloc2,js2))
     &                 * soangl(l,m2,js2,l,m1,js1)
                      sb(jsloc2)= sb(jsloc2) + 
     &                   CONJG(bcof2(band2,lm2,na,jsloc2,js2))
     &                 * soangl(l,m2,js2,l,m1,js1)

                    ENDDO ! m2  
                  ENDDO   ! jsloc2

                  DO jsloc1= 1,2
                    ra(jsloc1)= CMPLX(0.,0.) 
                    rb(jsloc1)= CMPLX(0.,0.) 
                    DO jsloc2= 1,2

                      ra(jsloc1)= ra(jsloc1) +  
     &                   sa(jsloc2) * rsopp(n,l,jsloc1,jsloc2) 
     &                 + sb(jsloc2) * rsoppd(n,l,jsloc1,jsloc2) 
                      rb(jsloc1)= rb(jsloc1) +
     &                   sa(jsloc2) * rsopdp(n,l,jsloc1,jsloc2)
     &                 + sb(jsloc2) * rsopdpd(n,l,jsloc1,jsloc2)
                       
                    ENDDO ! jsloc2
                  ENDDO   ! jsloc1 

                  DO jsloc1= 1,2
                    DO band1= 1,neigd

                      matel(band1,band2,n)= matel(band1,band2,n) +
     &                   acof1(band1,lm1,na,jsloc1,js1) * ra(jsloc1)   
     &                 + bcof1(band1,lm1,na,jsloc1,js1) * rb(jsloc1)   

                    ENDDO ! band1
                  ENDDO   ! jsloc1 

                ENDDO ! band2
              ENDDO   ! m1,lm1 
            ENDDO     ! l
          ENDDO       ! nat,na 
        ENDDO         ! n        
      ENDDO           ! js2,js1

      DO n= 1,ntype 
        IF ( (.not. soc_opt(ntype+1)) .or. soc_opt(n) ) THEN 
          DO band2= 1,neigd 
            DO band1= 1,neigd 
              matel(band1,band2,0)= matel(band1,band2,0) +
     &         matel(band1,band2,n) 
            ENDDO
          ENDDO 
        ENDIF 
      ENDDO 

      END SUBROUTINE ssomatel

