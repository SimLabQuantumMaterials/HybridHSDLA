      MODULE m_mpimakegroups
      CONTAINS
      SUBROUTINE mpi_make_groups(
     >                           irank,isize,nbasfcn,nkpt,neigd,
     >                           zelec,nw,nwd,ntypd,nvd,lmaxd,lmd,
     >                           nlod,llod,mlotot,mlolotot,jspd,n_u,
     >                           lmplmd,l_noco,l_ss,l_soc,l_constr,
     <                           n_start,n_groups,n,matsize,ne,
     <                           n_rank,n_size,SUB_COMM)
c------------------------------------------------------------------------
c
c Distribute the k-point / eigenvector  parallelisation so, that
c all pe's have aproximately equal load. Maximize for k-point 
c parallelisation. The naming conventions are as follows:
c
c groups             1               2          (n_groups = 2) 
c                 /     \         /     \
c k-points:      1       2       3       4      (nkpts = 4)
c               /|\     /|\     /|\     /|\    
c irank        0 1 2   3 4 5   0 1 2   3 4 5    (isize = 6)
c
c n_rank       0 1 2   0 1 2   0 1 2   0 1 2    (n_size = 3)
c
c nrec         1 2 3   4 5 6   7 8 9  10 11 12  ...rec. no. on eig-file
c              * *     * *     * *     *  *
c
c In the above example, 6 pe's should work on 4 k-points and distribute
c their load in a way, that 3 pe's work on each k-points, so 2 k-points
c are done in parellel (n_members=2) and there are 2 groups of k-points.
c n_rank and n_size are the equivalents of irank and isize. The former
c belong to the communicator SUB_COMM, the latter to MPI_COMM_WORLD.
c
c The results (containing a subset of ew & ev's) are written as separate
c records on the file 'eig' with the vacuum energy parameter of the
c marked (*) records set to 999.9 to indicate that this is only one part
c of a k-point's results. This marker is needed in routines fermie and 
c cdnval.
c          G.B. `99
c
c------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nbasfcn,nkpt,neigd,nw,nwd,irank,isize
      INTEGER, INTENT (IN)  :: ntypd,nvd,lmaxd,lmd,nlod,llod,mlotot
      INTEGER, INTENT (IN)  :: mlolotot,jspd,n_u,lmplmd
      LOGICAL, INTENT (IN)  :: l_noco,l_ss,l_soc,l_constr
      REAL,    INTENT (IN)  :: zelec
      INTEGER, INTENT (OUT) :: n_start,n_groups,n,SUB_COMM
      INTEGER, INTENT (OUT) :: matsize,n_rank,n_size,ne

      INTEGER i,n_members
      INTEGER, ALLOCATABLE :: i_mygroup(:)

      INCLUDE 'mpif.h'
      INTEGER WORLD_GROUP,SUB_GROUP
      INTEGER ierr(3)
      LOGICAL l_cm
c
      INTRINSIC nint, int, min
c     
c first determine the number of groups of k-points to process
c
      n_groups = 0
      IF (nkpt.GT.isize) THEN           ! if there are more k-points than PEs

        IF (mod(nkpt,isize).EQ.0) THEN  ! maybe nkpt is a multiple of isize
          n_groups = nkpt/isize
          n_size = 1
        ELSE                            ! or an integer fraction of isize fits
          DO i=2,isize
            IF (mod(isize,i).EQ.0) THEN
              n_size = i
              n_groups = nkpt * i/isize
              IF (mod(nkpt,isize/i).EQ.0) GOTO 990
            ENDIF
          ENDDO
          n_groups = nkpt               ! or use all PE's per k-point
          n_size = isize
        ENDIF

      ELSEIF (nkpt.LT.isize) THEN       ! if there are more PEs than k-points

        IF (mod(isize,nkpt).EQ.0) THEN  ! maybe isize is a multiple of nkpt 
           n_groups = 1
           n_size = isize/nkpt
        ELSE                            ! or an integer fraction of nkpts fits
          DO i=nkpt-1,2,-1
            IF (mod(nkpt,i).EQ.0) THEN
               n_groups = nkpt/i
               n_size = isize/i
               IF (mod(isize,i).EQ.0) GOTO 990
            ENDIF
          ENDDO
          n_groups = nkpt               ! or use all PE's per k-point
          n_size = isize
        ENDIF

      ELSE
!
! if there are as many pe's as k-points (isize = nkpt), use one PE per
! kpoint (n_size = 1) and finish in a single k-loop (n_groups = 1)
!
        n_groups = 1
        n_size = 1
      ENDIF



 990  IF (n_groups.EQ.0) STOP 'mpi_make_groups:1'
      n_members = nkpt/n_groups 
!
! check different algorithm
!
      CALL check_memory(
     >                  nbasfcn,neigd,ntypd,nvd,lmaxd,lmd,nlod,
     >                  llod,mlotot,mlolotot,jspd,n_u,lmplmd,
     >                  l_noco,l_ss,l_soc,l_constr,nkpt,isize,
     <                  n_size)

      write(*,*) n_size
      n_members = MIN(nkpt,isize)
      n_members = MIN(n_members , CEILING(REAL(isize)/n_size) ) + 1
      
      l_cm = .false.
      DO WHILE (.not.l_cm) 
        n_members = n_members - 1
        IF ((mod(isize,n_members) == 0).AND.
     +      (mod(nkpt,n_members) == 0) ) THEN
           l_cm = .true.
        ENDIF
      ENDDO
      n_groups = nkpt/n_members
      n_size = isize/n_members
      IF (irank == 0) THEN
        write(*,*) 'k-points in parallel: ',n_members
        write(*,*) "pe's per k-point:     ",n_size
        write(*,*) '# of k-point loops:   ',n_groups
      ENDIF
c
c now, we make the groups
c
      n_start = mod(irank,n_members) + 1
!!      n_start = INT(irank/n_size) * n_size
      ALLOCATE ( i_mygroup(n_size) )
      n = 0
      DO i = n_start,isize,n_members
!!      DO i = n_start+1,n_start+n_size
        n = n+1
        i_mygroup(n) = i-1
      ENDDO

c      write (*,*) irank,n_groups,n_start,i_mygroup

      CALL MPI_COMM_GROUP (MPI_COMM_WORLD,WORLD_GROUP,ierr)
      CALL MPI_GROUP_INCL (WORLD_GROUP,n_size,i_mygroup,
     +                                      SUB_GROUP,ierr)
      CALL MPI_COMM_CREATE (MPI_COMM_WORLD,SUB_GROUP,SUB_COMM,ierr) 

      CALL MPI_COMM_RANK (SUB_COMM,n_rank,ierr)
c
c determine number of columns per group
c
      n=0
      DO  i=1+n_rank, nbasfcn, n_size
        n=n+1
      ENDDO
      IF (n_size.GT.1) THEN
        matsize = nbasfcn * n
      ELSE
        matsize = (nbasfcn+1)*nbasfcn/2
      ENDIF
c
      IF (nw.EQ.nwd) THEN
        ne = neigd
      ELSE
        ne = nint( zelec/2. )
      ENDIF
      ne = max(ne,5)

      DEALLOCATE (i_mygroup)

      END SUBROUTINE mpi_make_groups

!----------------------------------------------------------------------
      SUBROUTINE check_memory(
     >                        nbasfcn,neigd,ntypd,nvd,lmaxd,lmd,nlod,
     >                        llod,mlotot,mlolotot,jspd,n_u,lmplmd,
     >                        l_noco,l_ss,l_soc,l_constr,nkpt,isize,
     <                        n_size)

!
! check the free and the (approximate) required memory ;
! determine minimal n_size to fit into the memory (hopefully).
!
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nbasfcn,neigd,ntypd,nvd,lmaxd,lmd,nlod
      INTEGER, INTENT (IN) :: llod,mlotot,mlolotot,jspd,n_u,lmplmd 
      INTEGER, INTENT (IN) :: nkpt,isize
      LOGICAL, INTENT (IN) :: l_noco,l_ss,l_soc,l_constr
      INTEGER, INTENT (OUT) :: n_size

      INTEGER  err, mb
      INTEGER*8 mem, matsz, m_h
      REAL, ALLOCATABLE :: test(:)

      n_size = CEILING( real(isize)/min(nkpt,isize) )
      
 10   CONTINUE
!
! some basic arrays allocated in eigen()
!
      mem = (lmplmd+1)*ntypd*4                       ! tuu,tdd etc.
      mem = mem + (lmd+1)*(2*llod+1)*max(mlotot,1)*2 ! tuulo ...
      mem = mem + (2*llod+1)**2 * max(mlolotot,1)    ! tuloulo
      IF (l_noco) mem = mem * 2                      ! both spins
      mem = mem + 49*n_u*jspd*2                      ! lda+U, *2 for complex
      mem = mem+INT((nbasfcn*2+(lmd+1)*ntypd)*0.5)+1 ! ind, *0.5 for integer

      matsz = nbasfcn * CEILING(REAL(nbasfcn)/n_size) ! size of a, b
#ifdef CPP_INVERSION
      mem = mem + 2 * matsz                          ! real case
#else
      mem = mem + 2 * matsz * 2                      ! complec case
#endif
!
! now the arrays in hssphn()
!
      m_h = nvd*(lmd+1)*4  + nvd*8 + nlod            ! ar, ai ..., cph, rph, vk, gk
      m_h = m_h + 2 * (2*llod+1)**2 * nlod * 3 * 2   ! alo,blo,clo
      IF (l_ss) m_h = m_h * 2
      m_h = m_h + nvd*(5+lmaxd)                      ! axr, ... plegend
      IF (l_ss.OR.l_constr.OR.(l_noco.AND.l_soc)) THEN
        m_h = m_h + nvd*(lmaxd+1)*ntypd*2*2          ! fj,gj
      ELSE 
        m_h = m_h + nvd*(lmaxd+1)*ntypd*2            ! fj,gj
      ENDIF
      IF (l_noco.AND.l_soc) THEN
        m_h = m_h + nvd*(lmaxd+4)
      ENDIF
      IF (l_constr) THEN
        m_h = m_h + (lmaxd+1)*ntypd
      ENDIF
      IF (l_noco.AND.(.NOT.l_ss)) THEN
        matsz = (nvd+mlotot) * CEILING(REAL(nvd+mlotot)/n_size)
        m_h = m_h + matsz * 2 * 2                    ! aahlp,bbhlp
      ENDIF
!
! see, whether it fits
!
      mb = (mem+m_h)*8/(1024)**2
      ALLOCATE ( test(mem+m_h) , stat = err)
      WRITE(*,*) mb,'Mbytes needed  in hssphn!',err,mem
      IF ( err /= 0 ) THEN
        n_size = n_size * 2
        IF (n_size > isize) THEN
          mb = (mem+m_h)*8/(1024)**2
          WRITE(*,*) mb,'Mbytes needed  in hssphn!'
          STOP 'mpi_make_groups: memory too small!'
        ENDIF 
        GOTO 10
      ENDIF
      DEALLOCATE (test)
!
! now, allocate z and jump into chani
!      
      matsz = nbasfcn * CEILING(REAL(nbasfcn)/n_size)   ! size of z
#ifdef CPP_INVERSION
      mem = mem + matsz                             ! real case
#else
      mem = mem + matsz * 2                         ! complex case
#endif
      mem = mem + matsz * 2 * 3                     ! asca,bsca,eigvec
      mem = mem + nbasfcn 
#ifdef CPP_INVERSION
      mem = mem + matsz                             ! work
#else
      mem = mem + matsz * 2                         ! work, rwork neglected
#endif
!
! see, whether it fits
!
      mb = (mem)*8/(1024)**2
      ALLOCATE ( test(mem) , stat = err)
      WRITE(*,*) mb,'Mbytes needed  in chani !',err,mem
      IF ( err /= 0 ) THEN
        n_size = n_size * 2
        IF (n_size > isize) THEN
          mb = (mem)*8/(1024)**2
          WRITE(*,*) mb,'Mbytes needed  in chani !'
          STOP 'mpi_make_groups: memory too small!'
        ENDIF
        GOTO 10
      ENDIF
      DEALLOCATE (test)
      
      END SUBROUTINE check_memory
!----------------------------------------------------------------------

      END MODULE m_mpimakegroups
