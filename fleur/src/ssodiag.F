PROGRAM ssodiag

IMPLICIT NONE

LOGICAL :: merging, matform, allk  

#ifdef CPP_MPI
INCLUDE 'mpif.h'
INTEGER :: ierr(3), stt(MPI_STATUS_SIZE)
#endif
INTEGER :: TYP_REAL, irank, isize, nnode, infile, outfile   
CHARACTER(len=9) :: matfile 

INTEGER :: matelversion 
REAL    :: timesum(2), time1, time2 
LOGICAL :: tkbin, eachtype, allocc, diagall    
INTEGER :: i1, i2, cyc, ncyc, ntypso, ntyp, cycnode, ncycnode 
INTEGER :: kp, cyclen, band1, band2, neig 
REAL    :: r1, r2, r3, qss(3), theta, phi, zelec,  &
           tkb, tkb0, fermi, fermiso, limdiag(2,2), eigsum, eigsumso  
COMPLEX :: c1 
REAL,    ALLOCATABLE :: time(:,:), timenode(:,:) 
INTEGER, ALLOCATABLE :: nbas(:,:), nbasnode(:,:)  
INTEGER, ALLOCATABLE :: knr(:,:) 
REAL,    ALLOCATABLE :: bkpt(:,:,:), bkptnode(:,:,:), wkpt(:) 
REAL,    ALLOCATABLE :: eig(:,:,:), eignode(:,:,:), eigso(:,:), eigsonode(:,:) 
COMPLEX, ALLOCATABLE :: mats(:,:,:), matn(:,:,:)

INTEGER :: knext 


  matform= .false. 
  ! matform=.true. => "matel*" formatted


  CALL cpu_time(time1)

  matelversion= 1 

#ifdef CPP_MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD,isize,ierr)
#ifdef CPP_T90
    CALL MPI_TYPE_CONTIGUOUS(1,MPI_REAL,TYP_REAL,ierr)
#else
    CALL MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_PRECISION,TYP_REAL,ierr)
#endif
  CALL MPI_TYPE_COMMIT(TYP_REAL,ierr)
#else
  irank= 0 
  isize= 1
  TYP_REAL= 0 
#endif

  infile=  11 
  outfile= 12  

  IF (matform) THEN
    INQUIRE(file='ssopart_form',exist=allk) ;  allk= (.not. allk) 
  ELSE
    INQUIRE(file='ssopart',exist=allk) ;  allk= (.not. allk) 
  ENDIF 
  IF (allk) THEN
    OPEN(infile,file='ssoinp',form='formatted',action='read')
  ELSE
    IF (matform) THEN
      OPEN(infile,file='ssopart_form',form='formatted',action='read') 
    ELSE 
      OPEN(infile,file='ssopart',form='unformatted',action='read') 
    ENDIF
  ENDIF 

  IF (irank==0) THEN 

    IF (allk) THEN 
      DO i1= 1,4
        READ(infile,fmt='(x)')  
      ENDDO
      READ(infile,fmt='(11x,l1)') diagall
      DO i1= 1,2 
        READ(infile,fmt='(x)')  
        READ(infile,fmt='(x)')  
        READ(infile,fmt='(11x,f14.10)') limdiag(i1,1)
        limdiag(i1,2)= limdiag(i1,1) 
      ENDDO 
      DO i1= 1,3 
        READ(infile,fmt='(x)')  
      ENDDO
      READ(infile,fmt='(11x,l1)') tkbin  
      READ(infile,fmt='(x)')  
      READ(infile,fmt='(11x,f14.10)') tkb 
    ELSE 
      IF (matform) THEN 
        READ(infile,fmt='(x)')
        DO i1= 1,2
          READ(infile,fmt='(2(e24.18,x))') limdiag(i1,1), limdiag(i1,2)
        ENDDO
        READ(infile,fmt='(l1,2(x,e24.18))') diagall, fermi, tkb  
        tkbin= .false. 
      ELSE
        READ(infile)
        DO i1= 1,2
          READ(infile) limdiag(i1,1), limdiag(i1,2)
        ENDDO
        READ(infile) diagall, fermi, tkb  
        tkbin= .false. 
      ENDIF 
    ENDIF 

    WRITE(*,fmt='(x)') 
    IF (diagall) THEN
      WRITE(*,fmt='(A)') 'all eigenstates used for exact diagonalization'
    ELSE
      WRITE(*,fmt='(A,e24.18)')  'exact diagonalization  in interval fermi+-',limdiag(1,1)
      WRITE(*,fmt='(A,e24.18)')  '2nd order perturbation in interval fermi+-',limdiag(2,1)
    ENDIF
    WRITE(*,fmt='(x)') 

  ENDIF ! (irank==0) 

  CLOSE(infile) 

  IF (matform) THEN 
    OPEN(infile,file='matelsso_form',form='formatted',action='read')
  ELSE
    OPEN(infile,file='matelsso',form='unformatted',action='read')
  ENDIF 

  IF (irank==0) THEN 

    IF (matform) THEN 
      READ(infile,fmt='(x)')
      READ(infile,fmt='(9x,i4)') i2 
      READ(infile,fmt='(4i8,x,l1)') ncyc, cyclen, neig, ntypso, eachtype 
      READ(infile,fmt='(3(e24.18,x))') qss(1), qss(2), qss(3)
      READ(infile,fmt='(2(e24.18,x))') theta, phi   
      READ(infile,fmt='(2(e24.18,x))') zelec, tkb0  
      READ(infile,fmt='(l1,i4)') merging, i1
    ELSE
      READ(infile)
      READ(infile) i2 
      READ(infile) ncyc, cyclen, neig, ntypso, eachtype 
      READ(infile) qss(1), qss(2), qss(3)
      READ(infile) theta, phi   
      READ(infile) zelec, tkb0  
      READ(infile) merging, i1
    ENDIF 
    IF (i2/=matelversion) THEN
      WRITE(*,fmt='(A,i4,x,A,i4)')  &
       'matel* has the format version #',i2,';  this program expects #',matelversion  
      STOP 'Stop:  old version of "matel*"-file-format'
    ENDIF 
    IF (eachtype) THEN
      ntyp= ntypso
    ELSE
      ntyp= 0 
    ENDIF 
    IF ( (.not. merging) .and. (i1/=isize) ) THEN
      STOP 'Stop:  number of processors and "matel"-files do not match'
    ENDIF
    IF (tkbin) THEN 
      tkb= tkb0 
    ENDIF 
    IF (tkb<0.) THEN
      WRITE(*,fmt='(A)') 'no temperature broadening'
    ELSE 
      WRITE(*,fmt='(A,e24.18)')  'temperature broadening:   T k_B =',tkb 
    ENDIF 
    WRITE(*,fmt='(x)') 

    ALLOCATE( eig(neig,cyclen,ncyc), bkpt(3,cyclen,ncyc), wkpt(ncyc), knr(cyclen,ncyc) ) 
    ALLOCATE( nbas(3,ncyc), time(2,ncyc), eigso(neig*cyclen,ncyc) )

    DO cyc= 1,ncyc
      IF (matform) THEN
        READ(infile,fmt='(9x,e24.18)') wkpt(cyc)
      ELSE
        READ(infile) wkpt(cyc)
      ENDIF 
      DO kp= 1,cyclen
        IF (matform) THEN 
          READ(infile,fmt='(5x,i8)') knr(kp,cyc)
          DO band1= 1,neig
            READ(infile,fmt='(9x,e24.18)') eig(band1,kp,cyc) 
          ENDDO
        ELSE
          READ(infile) knr(kp,cyc)
          DO band1= 1,neig
            READ(infile) eig(band1,kp,cyc) 
          ENDDO
        ENDIF 
        DO band1= 1,neig 
          eigso((kp-1)*neig+band1,cyc)= eig(band1,kp,cyc)   
           ! eig does not have appropiate format for fermie 
        ENDDO 
      ENDDO
    ENDDO 


    IF (merging) THEN 

      DO nnode= 0, isize-1
        IF (nnode<ncyc) THEN
          matfile(1:8)= 'matelXXX' 
          matfile(6:6)= ACHAR(IACHAR('0')+nnode/100 )
          matfile(7:7)= ACHAR(IACHAR('0')+MOD(nnode,100)/10 )
          matfile(8:8)= ACHAR(IACHAR('0')+MOD(nnode,10) )
          OPEN(outfile+nnode,file=matfile(1:8),form='unformatted',status='replace') 
        ENDIF
      ENDDO

      IF (matform) THEN
        READ(infile,fmt='(x)')
      ENDIF  

      DO cyc= 1,ncyc
        nnode= MOD(cyc-1,isize) 

        DO kp= 1,cyclen
          IF (matform) THEN 
            READ(infile,fmt='(6x,i8,3(x,e24.18))') i1, r1, r2, r3 
          ELSE
            READ(infile) i1, r1, r2, r3 
          ENDIF 
          WRITE(outfile+nnode) i1, r1, r2, r3
          DO i1= 0,ntyp 
            DO band2= 1,neig
              DO band1= 1,band2+neig
                IF (matform .and. merging) THEN
                  READ(infile,fmt='(2(e24.18,x))') r1, r2
                  c1= CMPLX(r1,r2) 
                ELSE
                  READ(infile) c1
                ENDIF 
                WRITE(outfile+nnode) c1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
 
      ENDDO

      DO nnode= 0, isize-1
        IF (nnode<ncyc) THEN 
          CLOSE(outfile+nnode) ! file='matelXXX' 
        ENDIF 
      ENDDO 

    ENDIF ! (merging) 


    IF (allk) THEN 
      CALL fermie(tkb,ncyc,neig*cyclen,eigso,REAL(cyclen)*zelec,wkpt, fermi,eigsum,allocc) 
      IF (allocc) THEN
        STOP 'Stop:  no unoccupied states at tkb=0'
      ENDIF 
    ENDIF 

  ENDIF ! (irank==0) 

  CLOSE(infile) ! file='matelsso' 

#ifdef CPP_MPI 
  CALL MPI_BCAST(ncyc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  CALL MPI_BCAST(cyclen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  CALL MPI_BCAST(neig,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  CALL MPI_BCAST(ntyp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  CALL MPI_BCAST(diagall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(fermi,1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(limdiag,4,TYP_REAL,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(merging,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif 

  IF ( irank < ncyc ) THEN

    ncycnode= ncyc/isize
    IF ( MOD(ncyc,isize) > irank ) THEN
      ncycnode= ncycnode + 1
    ENDIF

    ALLOCATE( nbasnode(3,ncycnode), timenode(2,ncycnode), eigsonode(neig*cyclen,ncycnode) )
    ALLOCATE( bkptnode(3,cyclen,ncycnode) )
    ALLOCATE( eignode(neig,cyclen,ncycnode) )
    ALLOCATE( mats(neig,neig,cyclen), matn(neig,neig,cyclen) )

    ! send eig 
    IF (irank==0) THEN 
      DO cyc= 1,ncyc
        nnode= MOD(cyc-1,isize)
        IF (nnode==0) THEN
          cycnode= (cyc-1)/isize+1
          DO kp= 1,cyclen
            DO band1= 1,neig 
              eignode(band1,kp,cycnode)= eig(band1,kp,cyc)  
            ENDDO
          ENDDO
#ifdef CPP_MPI
        ELSE
          CALL MPI_SEND(eig(1,1,cyc),neig*cyclen,TYP_REAL,nnode,0,MPI_COMM_WORLD,ierr)
#endif 
        ENDIF
      ENDDO
    ENDIF 
#ifdef CPP_MPI
    IF (irank/=0) THEN
      DO cycnode= 1,ncycnode
        CALL MPI_RECV(eignode(1,1,cycnode),neig*cyclen,TYP_REAL,0,0,MPI_COMM_WORLD,stt,ierr)
      ENDDO
    ENDIF
#endif

    matfile(1:8)= 'matelXXX'
    matfile(6:6)= ACHAR(IACHAR('0')+irank/100 )
    matfile(7:7)= ACHAR(IACHAR('0')+MOD(irank,100)/10 )
    matfile(8:8)= ACHAR(IACHAR('0')+MOD(irank,10) )
    OPEN(infile,file=matfile(1:8),form='unformatted',position='rewind')

    DO cycnode= 1,ncycnode 

      cyc= (cycnode-1)*isize+irank+1 

      WRITE(*,fmt='(2(A,i4,x))') 'cyc =',cyc,'/',ncyc 

      DO kp= 1, cyclen 
        READ(infile)  &
         i1, bkptnode(1,kp,cycnode), bkptnode(2,kp,cycnode), bkptnode(3,kp,cycnode)  
        DO band2= 1,neig
          DO band1= 1,band2
            READ(infile) mats(band1,band2,kp)
            mats(band2,band1,kp)= CONJG(mats(band1,band2,kp))
          ENDDO
        ENDDO
        IF (eachtype) THEN
          DO i1=1,ntyp ; DO band2=1,neig ; DO band1=1,band2 
            READ(infile)
          ENDDO ; ENDDO ; ENDDO
        ENDIF 
        DO band2= 1,neig
          DO band1= 1,neig 
            READ(infile) matn(band1,band2,kp)
          ENDDO
        ENDDO
        IF (eachtype) THEN
          DO i1=1,ntyp*neig*neig 
            READ(infile)
          ENDDO 
        ENDIF 
      ENDDO 

      ! enhance soc for test
      ! DO kp=1,cyclen ; DO band2=1,neig ; DO band1=1,neig
      !   mats(band1,band2,kp)= mats(band1,band2,kp)*3.  
      !   matn(band1,band2,kp)= matn(band1,band2,kp)*3.  
      ! ENDDO ; ENDDO ; ENDDO 

      CALL energies( (isize==1),  &  
       neig,cyclen,eignode(1,1,cycnode),mats,matn,diagall,fermi,limdiag,  &
       nbasnode(1,cycnode),timenode(1,cycnode),eigsonode(1,cycnode) )

    ENDDO ! cycnode= 1,ncycnode  

    DEALLOCATE( eignode, mats, matn )

    IF (merging) THEN 
      CLOSE(infile,status='delete') ! file= 'matelXXX'  
    ELSE
      CLOSE(infile) ! file= 'matelXXX'  
    ENDIF 

    ! send nbasnode,timenode,eigsonode,bkptnode
#ifdef CPP_MPI
    IF (irank/=0) THEN
      DO cycnode= 1,ncycnode 
        CALL MPI_SEND(nbasnode(1,cycnode),3,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(timenode(1,cycnode),2,TYP_REAL,0,0,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(eigsonode(1,cycnode),neig*cyclen,TYP_REAL,0,0,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(bkptnode(1,1,cycnode),3*cyclen,TYP_REAL,0,0,MPI_COMM_WORLD,ierr)
      ENDDO
    ENDIF 
#endif
    IF (irank==0) THEN
      DO cyc= 1,ncyc
        nnode= MOD(cyc-1,isize)
        IF (nnode==0) THEN
          cycnode= (cyc-1)/isize+1 
          DO band1= 1,3 
            nbas(band1,cyc)= nbasnode(band1,cycnode)
          ENDDO 
          DO band1= 1,2 
            time(band1,cyc)= timenode(band1,cycnode) 
          ENDDO 
          DO band1= 1,neig*cyclen 
            eigso(band1,cyc)= eigsonode(band1,cycnode)
          ENDDO
          DO kp= 1,cyclen 
            DO i1= 1,3
              bkpt(i1,kp,cyc)= bkptnode(i1,kp,cycnode)
            ENDDO
          ENDDO
#ifdef CPP_MPI
        ELSE 
          CALL MPI_RECV(nbas(1,cyc),3,MPI_INTEGER,nnode,0,MPI_COMM_WORLD,stt,ierr)
          CALL MPI_RECV(time(1,cyc),2,TYP_REAL,nnode,0,MPI_COMM_WORLD,stt,ierr)
          CALL MPI_RECV(eigso(1,cyc),neig*cyclen,TYP_REAL,nnode,0,MPI_COMM_WORLD,stt,ierr)
          CALL MPI_RECV(bkpt(1,1,cyc),3*cyclen,TYP_REAL,nnode,0,MPI_COMM_WORLD,stt,ierr)
#endif
        ENDIF 
      ENDDO ! cyc= 1,ncyc  
    ENDIF 

    DEALLOCATE( nbasnode, timenode, eigsonode, bkptnode )

  ENDIF ! ( irank < ncyc )


  IF (irank==0) THEN 

    ! CALL bandplot(outfile,ncyc,cyclen,neig,neig*cyclen,bkpt,eig,eigso)  ! test  

    IF (allk) THEN 

      CALL fermie(tkb,ncyc,neig*cyclen,eigso,REAL(cyclen)*zelec,wkpt, fermiso,eigsumso,allocc)  
      IF (allocc) THEN
        STOP 'Stop:  no unoccupied states at tkb=0'
      ENDIF 

      DO i1= 1,2 
        timesum(i1)= 0. 
        DO cyc= 1,ncyc
          timesum(i1)= timesum(i1) + time(i1,cyc)
        ENDDO
      ENDDO 

      CALL cpu_time(time2) 

      OPEN(outfile,file='eigsso',form='formatted',status='replace')
      WRITE(outfile,fmt='(A,x,e24.18)')  'total sum of eigenvalues with SOC    =', eigsumso 
      WRITE(outfile,fmt='(A,x,e24.18)')  'total sum of eigenvalues without SOC =', eigsum   
      WRITE(outfile,fmt='(x)') 
      WRITE(outfile,fmt='(A,x,e24.18)')  'fermi-energy with SOC    =', fermiso 
      WRITE(outfile,fmt='(A,x,e24.18)')  'fermi-energy without SOC =', fermi 
      IF (diagall) THEN
        WRITE(outfile,fmt='(A)') 'all eigenstates used for exact diagonalization'
      ELSE 
        WRITE(outfile,fmt='(A,2(e24.18,A))')  'exact diagonalization  in interval [',  &
         fermi-limdiag(1,1),',',fermi+limdiag(1,2),']'
        WRITE(outfile,fmt='(A,2(e24.18,A))')  '2nd order perturbation in interval [',  &
         fermi-limdiag(2,1),',',fermi+limdiag(2,2),']'
      ENDIF 
      WRITE(outfile,fmt='(x)') 
      WRITE(outfile,fmt='(3(A,x,e24.18),x,A)')  'qss = (',qss(1),',',qss(2),',',qss(3),')'
      WRITE(outfile,fmt='(2(A,f6.3,x))') 'theta/pi =',theta,',  phi/pi =', phi
      WRITE(outfile,fmt='(x)') 
      IF (tkb<0.) THEN
        WRITE(outfile,fmt='(A)') 'no temperature broadening'
      ELSE
        WRITE(outfile,fmt='(A,e24.18)')  'temperature broadening:   T k_B =',tkb
      ENDIF 
      WRITE(outfile,fmt='(x)') 
      WRITE(outfile,fmt='(x)') 
      i1= INT(timesum(1)/10.)
      WRITE(outfile,fmt='(A,f8.1,x,A,i2,A,2i1,A,i1,f3.1,x,A)')  &
       'time for exact diag.  =',timesum(1),'s =',  &
       i1/360,':',MOD(i1,360)/60,MOD(i1,60)/6,':',MOD(i1,6),timesum(1)-REAL(10*i1),'h:min:s'
      i1= INT(timesum(2)/10.)
      WRITE(outfile,fmt='(A,f8.1,x,A,i2,A,2i1,A,i1,f3.1,x,A)')  &
       'time for perturbation =',timesum(2),'s =',  &
       i1/360,':',MOD(i1,360)/60,MOD(i1,60)/6,':',MOD(i1,6),timesum(2)-REAL(10*i1),'h:min:s'
      i1= INT((time2-time1)/10.)
      WRITE(outfile,fmt='(A,f8.1,x,A,i2,A,2i1,A,i1,f3.1,x,A,i3,x,A)')  &
       'time for everything   =',time2-time1,'s =',  &
       i1/360,':',MOD(i1,360)/60,MOD(i1,60)/6,':',MOD(i1,6),time2-time1-REAL(10*i1),  &
       'h:min:s   on each of',isize,'processors'
      WRITE(outfile,fmt='(x)') 
      WRITE(outfile,fmt='(x)') 
      DO cyc= 1,ncyc
        WRITE(outfile,fmt='(A,i4,2x,A,i5,2(A,i6),2(A,f8.2),A)') 'cyc =',cyc,  &
         '( nbas=',nbas(1,cyc),', nocc=',nbas(2,cyc),', nuno=',nbas(3,cyc),  &
         ',  dtime=',time(1,cyc),'s, ptime=',time(2,cyc),'s ) :' 
        WRITE(outfile,fmt='(x)')
        WRITE(outfile,fmt='(A,x,e24.18)') 'weight =',wkpt(cyc)
        WRITE(outfile,fmt='(A)') 'eigenvalues with SOC:'
        WRITE(outfile,fmt='(x,5f11.6)') (eigso(band1,cyc),band1=1,neig*cyclen)
        WRITE(outfile,fmt='(A)') 'eigenvalues without SOC:'
        DO kp= 1,cyclen
          WRITE(outfile,fmt='(i3,A,i4,x,A,3(e24.18,A))')  &
           kp,'. k-point #',knr(kp,cyc),'at (',  &
           bkpt(1,kp,cyc),',',bkpt(2,kp,cyc),',',bkpt(3,kp,cyc),')' 
          WRITE(outfile,fmt='(x,5f11.6)') (eig(band1,kp,cyc),band1=1,neig)
        ENDDO
        WRITE(outfile,fmt='(x)')
      ENDDO ! cyc= 1,ncyc
      CLOSE(outfile) 

      WRITE(*,fmt='(x)') 
      WRITE(*,fmt='(A,x,e24.18)')  'total sum of eigenvalues with SOC    =', eigsumso 
      WRITE(*,fmt='(A,x,e24.18)')  'total sum of eigenvalues without SOC =', eigsum   
      WRITE(*,fmt='(x)') 
      WRITE(*,fmt='(A,x,e24.18)')  'fermi-energy with SOC    =', fermiso 
      WRITE(*,fmt='(A,x,e24.18)')  'fermi-energy without SOC =', fermi 
      IF (diagall) THEN
        WRITE(*,fmt='(A)') 'all eigenstates used for exact diagonalization'
      ELSE 
        WRITE(*,fmt='(A,2(e24.18,A))')  'exact diagonalization  in interval [',  &
         fermi-limdiag(1,1),',',fermi+limdiag(1,2),']'
        WRITE(*,fmt='(A,2(e24.18,A))')  '2nd order perturbation in interval [',  &
         fermi-limdiag(2,1),',',fermi+limdiag(2,2),']'
      ENDIF 
      WRITE(*,fmt='(x)') 

    ELSE ! .not. (allk)  

      OPEN(outfile,file='ssopartso',form='unformatted',status='replace')
      WRITE(outfile) 
      WRITE(outfile) ncyc, cyclen, neig 
      WRITE(outfile) fermi, diagall, limdiag, tkb, zelec 
      WRITE(outfile) qss, theta, phi  
      DO cyc= 1,ncyc 
        WRITE(outfile) nbas(1,cyc),nbas(2,cyc),nbas(3,cyc)  
        WRITE(outfile) time(1,cyc),time(2,cyc) 
        WRITE(outfile) wkpt(cyc) 
        DO kp= 1,cyclen
          WRITE(outfile) bkpt(1,kp,cyc),bkpt(2,kp,cyc),bkpt(3,kp,cyc)
          DO band1= 1,neig
            WRITE(outfile) eig(band1,kp,cyc) 
            WRITE(outfile) eigso(band1+(kp-1)*neig,cyc) 
          ENDDO
        ENDDO
      ENDDO
      CLOSE(outfile) 

    ENDIF ! (allk) ELSE  

    DEALLOCATE( eig,bkpt,wkpt,knr, nbas,time,eigso )

  ENDIF ! (irank==0) 

#ifdef CPP_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF (irank==0) THEN
#endif
    WRITE(*,fmt='(A)') 'all done'
#ifdef CPP_MPI
  ENDIF 
  CALL MPI_TYPE_FREE(TYP_REAL,ierr) 
  CALL MPI_FINALIZE(ierr)
#endif

END PROGRAM

! ==============================================================================

SUBROUTINE bandplot(outfile,ncyc,cyclen,neig,ne,bkpt,eig,eigso)

IMPLICIT NONE 

INTEGER, INTENT(IN) :: outfile, ncyc, cyclen, neig
INTEGER, INTENT(IN) :: ne(ncyc)
REAL,    INTENT(IN) :: bkpt(3,cyclen,ncyc), eig(neig,cyclen,ncyc), eigso(neig*cyclen,ncyc)

INTEGER :: cyc, kp, band1
CHARACTER(len=19) :: fmttxt

  IF (neig>9999) THEN
    STOP 'Stop in bandplot:  neig*cyclen > 9999'
  ENDIF
  fmttxt='(e12.6,x,????f10.6)'
  fmttxt(10:10)= ACHAR(IACHAR('0')+(neig*cyclen)/1000 )
  fmttxt(11:11)= ACHAR(IACHAR('0')+MOD(neig*cyclen,1000)/100 )
  fmttxt(12:12)= ACHAR(IACHAR('0')+MOD(neig*cyclen,100)/10 )
  fmttxt(13:13)= ACHAR(IACHAR('0')+MOD(neig*cyclen,10) )
  OPEN(outfile,file='ebands',form='formatted',status='replace')
  DO cyc= 1,ncyc
    WRITE(outfile,fmt=fmttxt) bkpt(1,1,cyc), ((eig(band1,kp,cyc),band1=1,neig),kp=1,cyclen) 
  ENDDO 
  CLOSE(outfile)

  DO cyc= 1,ncyc
    IF (ne(cyc)/=ne(1)) THEN
      STOP 'Stop in bandplot:  ne not constant'
    ENDIF
  ENDDO

  IF (ne(1)>9999) THEN
    STOP 'Stop in bandplot:  ne(1) > 9999'
  ENDIF
  fmttxt='(e12.6,x,????f10.6)'
  fmttxt(10:10)= ACHAR(IACHAR('0')+ne(1)/1000 )
  fmttxt(11:11)= ACHAR(IACHAR('0')+MOD(ne(1),1000)/100 )
  fmttxt(12:12)= ACHAR(IACHAR('0')+MOD(ne(1),100)/10 )
  fmttxt(13:13)= ACHAR(IACHAR('0')+MOD(ne(1),10) )
  OPEN(outfile,file='esocbands',form='formatted',status='replace')
  DO cyc= 1,ncyc
    WRITE(outfile,fmt=fmttxt) bkpt(1,1,cyc), (eigso(band1,cyc),band1=1,ne(1))
  ENDDO 
  CLOSE(outfile)

END SUBROUTINE bandplot 

! ==============================================================================

SUBROUTINE energies(  &
 writing,neig,cyclen,eig,mats,matn,diagall,fermi,limdiag, nbasout,time,eigso)

IMPLICIT NONE

LOGICAL, INTENT(IN)  :: writing, diagall 
INTEGER, INTENT(IN)  :: neig, cyclen
REAL,    INTENT(IN)  :: eig(neig,cyclen), fermi, limdiag(2,2) 
COMPLEX, INTENT(IN)  :: mats(neig,neig,cyclen), matn(neig,neig,cyclen)

INTEGER, INTENT(OUT) :: nbasout(3)
REAL,    INTENT(OUT) :: time(2) 
REAL,    INTENT(OUT) :: eigso(neig*cyclen)

REAL    :: time1
INTEGER :: i1, kp, km, band1, band2, nbas, nocc, nuno 
INTEGER :: bidx(neig*cyclen,-1:1), kidx(neig*cyclen,-1:1)
REAL    :: r1
COMPLEX :: c1, c2  
REAL,    ALLOCATABLE :: eigdiag(:) 
COMPLEX, ALLOCATABLE :: matbig(:,:)

  DO i1= 1,2
    IF (.not. diagall) THEN
      IF (limdiag(2,i1)+1.e-10 > limdiag(1,i1)) THEN
        WRITE(*,fmt='(x)') 
        WRITE(*,fmt='(2(A,i1,A,x,e24.18,2x))')  &
         'WARNING:   limdiag(1,',i1,') =',limdiag(1,i1),'limdiag(2,',i1,') =',limdiag(2,i1) 
        WRITE(*,fmt='(x)') 
      ENDIF
    ENDIF 
  ENDDO 

  nbas= 0
  nocc= 0
  nuno= 0  
  DO kp= 1,cyclen
    DO band1= 1,neig
      IF ( diagall .or.                                      &
           ( ( eig(band1,kp) > fermi-limdiag(1,1) ) .and.    & 
             ( eig(band1,kp) < fermi+limdiag(1,2) )       )  &
         ) THEN  
        nbas= nbas+1
        bidx(nbas,0)= band1 
        kidx(nbas,0)= kp  
      ELSE
        IF ( eig(band1,kp) < fermi ) THEN  
          nocc= nocc+1
          bidx(nocc,-1)= band1
          kidx(nocc,-1)= kp
        ELSE
          nuno= nuno+1
          bidx(nuno,+1)= band1
          kidx(nuno,+1)= kp
        ENDIF  
      ENDIF  
    ENDDO
  ENDDO 

  nbasout(1)= nbas
  nbasout(2)= nocc
  nbasout(3)= nuno 
  IF (writing) THEN 
    WRITE(*,fmt='(A,i5,2(x,A,i6))') 'nbas =',nbas,',  nocc =',nocc,',  nuno =',nuno 
  ENDIF 
  time(1)= 0.  ;  time(2)= 0. 

  IF (nbas>0) THEN 

    ALLOCATE( matbig(nbas,nbas), eigdiag(nbas) )

    DO band1= 1,nbas

      DO band2= 1,band1 
        matbig(band1,band2)= CMPLX(0.,0.)
        IF ( ( ABS(kidx(band1,0)-kidx(band2,0)) <= 1 ) .or.          &
             ( ABS(kidx(band1,0)-kidx(band2,0)) == cyclen-1 ) ) THEN
          IF ( kidx(band1,0)==kidx(band2,0) ) THEN
            matbig(band1,band2)= mats(bidx(band1,0),bidx(band2,0),kidx(band1,0))
          ENDIF 
          IF ( MOD(kidx(band1,0),cyclen)+1==kidx(band2,0) ) THEN
            matbig(band1,band2)= matbig(band1,band2) +  & 
             matn(bidx(band1,0),bidx(band2,0),kidx(band1,0)) 
          ENDIF 
          IF ( kidx(band1,0)==MOD(kidx(band2,0),cyclen)+1 ) THEN
            matbig(band1,band2)= matbig(band1,band2) +  &
             CONJG(matn(bidx(band2,0),bidx(band1,0),kidx(band2,0)))  
          ENDIF
        ENDIF
        matbig(band2,band1)= CONJG(matbig(band1,band2))
      ENDDO

      matbig(band1,band1)=  &
       matbig(band1,band1) + CMPLX(eig(bidx(band1,0),kidx(band1,0)),0.)

    ENDDO

    IF (writing) THEN 
      WRITE(*,fmt='(A)') 'Diagonalizing ...'
    ENDIF 
    CALL cpu_time(time1)
    CALL eigvas(nbas,diagall, matbig, eigdiag)
    CALL cpu_time(time(1)) 
    time(1)=time(1)-time1 
    IF (writing) THEN
      WRITE(*,fmt='(A)') ' ... done'
    ENDIF 

    DO band1= 1,nbas
      eigso(nocc+band1)= eigdiag(band1) 
    ENDDO 

  ENDIF ! (nbas>0) 

  IF (nocc+nuno > 0) THEN

    IF (nocc > 0) THEN 
      DO band1= 1,nocc
        eigso(band1)= eig(bidx(band1,-1),kidx(band1,-1))
      ENDDO
    ENDIF 
    IF (nuno > 0) THEN
      DO band1= 1,nuno
        eigso(nocc+nbas+band1)= eig(bidx(band1,+1),kidx(band1,+1))
      ENDDO
    ENDIF 

    ! GOTO 100 ! test 

    IF (writing) THEN 
      WRITE(*,fmt='(A)') 'RS-perturbation ...'
    ENDIF 
    CALL cpu_time(time1) 

    ! 1st order pert
    IF (nocc > 0) THEN 
      DO band1= 1,nocc
        eigso(band1)= eigso(band1) +  &
         REAL(mats(bidx(band1,-1),bidx(band1,-1),kidx(band1,-1)))  
        IF (cyclen==1) THEN
          eigso(band1)= eigso(band1) +  &
           2.*REAL(matn(bidx(band1,-1),bidx(band1,-1),kidx(band1,-1)))  
        ENDIF 
      ENDDO
    ENDIF

    ! 2nd order pert 
    IF ( (nbas>0) .and. (nocc>0) ) THEN
      DO band2= 1,nbas
        IF ( eigdiag(band2) > fermi-limdiag(2,1) ) THEN 
          DO band1= 1,nocc
            c2= CMPLX(0.,0.) 
            DO i1= 1,nbas
              IF ( ( ABS(kidx(i1,0)-kidx(band1,-1)) <= 1 ) .or.          &
                   ( ABS(kidx(i1,0)-kidx(band1,-1)) == cyclen-1 ) ) THEN
                c1= CMPLX(0.,0.) 
                IF ( kidx(i1,0)==kidx(band1,-1) ) THEN
                  c1= c1 + mats(bidx(i1,0),bidx(band1,-1),kidx(i1,0)) 
                ENDIF
                IF ( MOD(kidx(i1,0),cyclen)+1==kidx(band1,-1) ) THEN
                  c1= c1 + matn(bidx(i1,0),bidx(band1,-1),kidx(i1,0)) 
                ENDIF
                IF ( kidx(i1,0)==MOD(kidx(band1,-1),cyclen)+1 ) THEN
                  c1= c1 + CONJG(matn(bidx(band1,-1),bidx(i1,0),kidx(band1,-1))) 
                ENDIF
                c2= c2 + c1*CONJG(matbig(i1,band2))
              ENDIF 
            ENDDO
            r1=  REAL( CONJG(c2)*c2 )                                  &  
                /( eigdiag(band2)-eig(bidx(band1,-1),kidx(band1,-1)) )  
            eigso(nocc+band2)= eigso(nocc+band2) + r1 
            eigso(band1)= eigso(band1) - r1 
          ENDDO
        ENDIF
      ENDDO
    ENDIF ! (nbas>0) .and. (nocc>0) 
    IF ( (nbas>0) .and. (nuno>0) ) THEN
      DO band2= 1,nbas
        IF ( eigdiag(band2) < fermi+limdiag(2,2) ) THEN
          DO band1= 1,nuno
            c2= CMPLX(0.,0.)
            DO i1= 1,nbas
              IF ( ( ABS(kidx(i1,0)-kidx(band1,+1)) <= 1 ) .or.          &
                   ( ABS(kidx(i1,0)-kidx(band1,+1)) == cyclen-1 ) ) THEN
                c1= CMPLX(0.,0.)
                IF ( kidx(i1,0)==kidx(band1,+1) ) THEN
                  c1= c1 + mats(bidx(i1,0),bidx(band1,+1),kidx(i1,0))
                ENDIF
                IF ( MOD(kidx(i1,0),cyclen)+1==kidx(band1,+1) ) THEN
                  c1= c1 + matn(bidx(i1,0),bidx(band1,+1),kidx(i1,0))
                ENDIF
                IF ( kidx(i1,0)==MOD(kidx(band1,+1),cyclen)+1 ) THEN
                  c1= c1 + CONJG(matn(bidx(band1,+1),bidx(i1,0),kidx(band1,+1)))
                ENDIF
                c2= c2 + c1*CONJG(matbig(i1,band2))
              ENDIF
            ENDDO
            r1=  REAL( CONJG(c2)*c2 )                                  &                         
                /( eigdiag(band2)-eig(bidx(band1,+1),kidx(band1,+1)) )
            eigso(nocc+band2)= eigso(nocc+band2) + r1      
          ENDDO
        ENDIF
      ENDDO
    ENDIF ! (nbas>0) .and. (nuno>0) 
    IF ( (nocc>0) .and. (nuno>0) ) THEN
      DO band2= 1,nocc
        DO band1= 1,nuno
          c2= CMPLX(0.,0.)
          IF ( ( ABS(kidx(band2,-1)-kidx(band1,+1)) <= 1 ) .or.          &
               ( ABS(kidx(band2,-1)-kidx(band1,+1)) == cyclen-1 ) ) THEN
            IF ( kidx(band2,-1)==kidx(band1,+1) ) THEN
              c2= c2 + mats(bidx(band2,-1),bidx(band1,+1),kidx(band2,-1))
            ENDIF
            IF ( MOD(kidx(band2,-1),cyclen)+1==kidx(band1,+1) ) THEN
              c2= c2 + matn(bidx(band2,-1),bidx(band1,+1),kidx(band2,-1))
            ENDIF
            IF ( kidx(band2,-1)==MOD(kidx(band1,+1),cyclen)+1 ) THEN
              c2= c2 + CONJG(matn(bidx(band1,+1),bidx(band2,-1),kidx(band1,+1)))
            ENDIF
          ENDIF
          r1=  REAL( CONJG(c2)*c2 )                    &
              /(  eig(bidx(band2,-1),kidx(band2,-1))   & 
                 -eig(bidx(band1,+1),kidx(band1,+1)) )
          eigso(band2)= eigso(band2) + r1
        ENDDO
      ENDDO
    ENDIF ! (nocc>0) .and. (nuno>0)

    CALL cpu_time(time(2))
    time(2)= time(2)-time1 
    IF (writing) THEN 
      WRITE(*,fmt='(A)') ' ... done'
    ENDIF 

    100 CONTINUE 

  ENDIF ! (nocc+nuno > 0) 

  IF (nbas>0) THEN
    DEALLOCATE( matbig, eigdiag )
  ENDIF 

END SUBROUTINE energies

! ==============================================================================

SUBROUTINE eigvas(matsize,diagall, matrix, eigvals)

IMPLICIT NONE

LOGICAL, INTENT(IN)    :: diagall 
INTEGER, INTENT(IN)    :: matsize
REAL,    INTENT(OUT)   :: eigvals(matsize)
COMPLEX, INTENT(INOUT) :: matrix(matsize,matsize)

CHARACTER :: eigversion 
INTEGER   :: i1, i2
CHARACTER :: vectors
INTEGER   :: info
INTEGER   :: lrwork, lcwork, liwork
INTEGER, ALLOCATABLE :: isuppz(:)
COMPLEX, ALLOCATABLE :: matout(:,:)
COMPLEX, ALLOCATABLE :: cwork(:)
REAL,    ALLOCATABLE :: rwork(:)
INTEGER, ALLOCATABLE :: iwork(:)

#ifdef CPP_DOUBLE
# define lapack_zheev zheev 
# define lapack_zheevr zheevr 
#else
# define lapack_zheev cheev 
# define lapack_zheevr cheevr 
#endif 

  eigversion= '0'  ! zheev
  ! eigversion= 'R'  ! zheevr

  IF (diagall) THEN 
    vectors= 'N'
  ELSE
    vectors= 'V'
  ENDIF 

  IF (eigversion=='R') THEN

#ifdef CPP_DEC
    eigversion= '-' 
#else

    ALLOCATE( matout(matsize,matsize), isuppz(2*matsize) )

    lcwork= 1 ; lrwork= 1 ; liwork= 1

    DO i1= 1,2

      ALLOCATE( cwork(lcwork), rwork(lrwork), iwork(liwork) )

      IF (i1==1) THEN
        lcwork= -1  ;  lrwork= -1  ;  liwork= -1
      ENDIF

      CALL lapack_zheevr(  &
       vectors,'A','U',matsize,  &
       matrix,matsize,  &
       0.,1.,1,1,-1.,i2,  &
       eigvals,matout,matsize,  &
       isuppz,cwork,lcwork,rwork,lrwork,iwork,liwork,  &
       info )
      IF (info/=0) THEN
        STOP 'Stop in eigvas:  error while diagonalizing'
      ENDIF

      IF (i1==1) THEN
        lcwork= NINT(REAL(cwork(1))) ; lrwork= NINT(rwork(1)) ; liwork= iwork(1)
      ENDIF
  
      DEALLOCATE( cwork, rwork, iwork )

    ENDDO

    DO i2=1,matsize
      DO i1=1,matsize
        matrix(i1,i2)= matout(i1,i2)
      ENDDO
    ENDDO

    DEALLOCATE( matout, isuppz )

#endif 

  ENDIF     

  IF ( (eigversion=='0') .or. (eigversion=='-') ) THEN 

    IF (eigversion=='-') THEN
      WRITE (*,fmt='(x)') 
      WRITE (*,fmt='(A)') 'WARNING in eigvas:'  
      WRITE (*,fmt='(A)')  &
       'designated Lapack-routine not installed, using zheev instead!' 
      WRITE (*,fmt='(x)') 
    ENDIF 

    ALLOCATE( rwork(3*matsize), cwork(2*matsize) ) 

    CALL lapack_zheev(                             &
     vectors,'U',matsize,                          &
     matrix,                                       & 
     matsize,eigvals,cwork,2*matsize,rwork,info )
    IF (info/=0) THEN
      STOP 'Stop in eigvas:  error while diagonalizing'
    ENDIF

    DEALLOCATE( rwork, cwork )

  ENDIF 

END SUBROUTINE eigvas 

! ==============================================================================

SUBROUTINE fermie (tkb,nkpt,neigd,eig,zelec,wkpt, fermi,eigsum,allocc)

IMPLICIT NONE

REAL,    INTENT(IN)  :: tkb
INTEGER, INTENT(IN)  :: nkpt, neigd
REAL,    INTENT(IN)  :: eig(neigd,nkpt)  
REAL,    INTENT(IN)  :: zelec, wkpt(nkpt)

REAL,    INTENT(OUT) :: fermi, eigsum  
LOGICAL, INTENT(OUT) :: allocc 

INTEGER :: kp, je, nall, steps 
REAL    :: wsum, wone, wkmin, eps0, eps, ff, sff, dff, sdff   
REAL,    ALLOCATABLE :: eigall(:), wkall(:) 
INTEGER, ALLOCATABLE :: ne(:)

REAL  fermfun 

  eps0= 1.e-10 
  eps= 1.e-15 

  ALLOCATE( ne(nkpt) ) 
  DO kp= 1,nkpt
    ne(kp)= neigd
  ENDDO 

  nall= 0 
  DO kp= 1,nkpt
    nall= nall+ne(kp)
  ENDDO
  ALLOCATE( eigall(nall), wkall(nall) )

  wsum= 0. 
  DO kp=1,nkpt
    wsum= wsum + zelec*wkpt(kp)
  ENDDO 

  nall= 0
  DO kp= 1,nkpt
    IF (ne(kp)>0) THEN
      DO je= 1,ne(kp)
        nall= nall+1
        eigall(nall)= eig(je,kp)
        wkall(nall)= wkpt(kp)/wsum
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE( ne ) 

  CALL sort(nall,eigall,wkall) 

  wkmin= wkall(1)
  DO je= 1,nall
    IF ( wkmin > wkall(je) ) THEN
      wkmin= wkall(je)
    ENDIF
  ENDDO 
     
  je= 0  
  wone= 0. 
  eigsum= 0. 
  DO WHILE ( (wone+wkmin/2. < 1.) .and. (je < nall) ) 
    je= je+1
    wone= wone + wkall(je) 
    eigsum= eigsum + wkall(je)*eigall(je) 
  ENDDO 

  IF (je==nall) THEN

    allocc= .true.
    fermi= eigall(je)

  ELSE

    allocc= .false. 
    fermi= (eigall(je)+eigall(je+1))/2. 

    IF (tkb > 0.) THEN

      sff= 2. 
      steps= 0 
      DO WHILE ( ABS(sff-1.) > eps0 )
        steps= steps+1
        IF (steps>50) THEN
          STOP 'Stop:  determination of fermi energy not converged'
        ENDIF

        sff= 0.
        sdff= 0.
        DO je= 1,nall
          ff= fermfun(eigall(je),fermi,tkb) 
          dff= ff*(1.-ff)/tkb 
          sff= sff + ff*wkall(je)
          sdff= sdff + dff*wkall(je)
        ENDDO
        IF (abs(sdff)<eps) THEN 
          STOP 'Stop:  minimum of fermi energy to shallow'
        ENDIF
      
        fermi= fermi - (sff-1.)/sdff  

      ENDDO

      eigsum= 0. 
      DO je= 1,nall
        eigsum= eigsum + wkall(je)*eigall(je)*fermfun(eigall(je),fermi,tkb)
      ENDDO  

    ENDIF 

  ENDIF 

  eigsum= eigsum*wsum  

  DEALLOCATE( eigall, wkall ) 

END SUBROUTINE fermie

! ==============================================================================

SUBROUTINE sort (n, skio,weio)

IMPLICIT NONE

INTEGER, INTENT(IN)    :: n

REAL,    INTENT(INOUT) :: skio(n), weio(n)  

INTEGER :: jx(n)
REAL    :: sk(n), we(n)  
REAL    :: eps, q
INTEGER :: i, ind, ir, j, l

DATA  eps/1.e-10/

  DO i= 1,n
    jx(i)= i
    sk(i)= skio(i)
    we(i)= weio(i)
  ENDDO

  l= n/2+1
  ir= n
  DO WHILE (ir>1)

    IF (l>1) THEN
      l= l-1
      ind= jx(l)
      q= sk(ind)
    ELSE
      ind= jx(ir)
      q= sk(ind)
      jx(ir)= jx(1)
      ir= ir-1
    ENDIF

    IF (ir==1) THEN

      jx(1)= ind

    ELSE

      i= l
      j= l+l

      DO WHILE (j<=ir) 
        
        IF (j<ir) THEN
          IF ( sk(jx(j+1)) - sk(jx(j)) > eps ) THEN
            j= j+1
          ENDIF
        ENDIF
        IF ( sk(jx(j)) - q > eps ) THEN
          jx(i)= jx(j)
          i= j
          j= j+j
        ELSE
          j= ir+1
        ENDIF

      ENDDO

      jx(i)= ind  
 
    ENDIF

  ENDDO 

  DO i=1,n
    skio(i)= sk(jx(i))
    weio(i)= we(jx(i)) 
  ENDDO

END SUBROUTINE sort 

! ==============================================================================

REAL FUNCTION fermfun (eig,fermi,tkb)

IMPLICIT NONE

REAL, INTENT(IN) :: eig, fermi, tkb

REAL interv

  interv= 8.*tkb

  IF ( (ABS(eig-fermi) < interv) .and. (tkb > 0.) ) THEN
    fermfun= 1./( EXP( (eig-fermi) / tkb ) + 1. )  
  ELSE
    IF ( eig > fermi ) THEN
      fermfun= 0.
    ELSE
      fermfun= 1.
    ENDIF
  ENDIF

END FUNCTION fermfun 

! ==============================================================================

SUBROUTINE cpu_time(time)

!**********************************************************************
!     this soubroutine determines the CPU-time in seconds the
!     program needed since it has been started. the result is
!     stored in the variable time.
!     on IBM-AIX using xl-fortran this is done with the mclock
!     command, which has an accuracy of 1/100 of a second.
!     on other machines this command probably needs to be replaced.
!                                             p.kurz  8.2.96
!**********************************************************************

      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL time
!     ..
!     .. External Functions ..
#if (defined(CPP_AIX)||defined(CPP_APC))
!---> IBM AIX or ATHLON PC
      INTEGER mclock
      EXTERNAL mclock
#elif defined (CPP_T90)
!---> Cray T90, J90
      REAL second
      EXTERNAL second
#elif defined (CPP_VPP)
!---> vpp300
      REAL Y1
      EXTERNAL clock
#else
!---> DEC ALPHA
      REAL etime
      EXTERNAL etime
#endif
!     ..

#ifdef CPP_AIX
!---> IBM RS6000
      time = mclock()*0.01
#elif defined (CPP_APC)
      time = mclock()*0.001
#elif defined (CPP_T90)
!---> Cray T90, J90
      time = second()
#elif defined (CPP_VPP)
!---> vpp300
      CALL CLOCK(Y1,2,2)
      time = Y1*0.000001
#else
!---> DEC ALPHA
      time = etime()
#endif

END SUBROUTINE cpu_time

! ==============================================================================

