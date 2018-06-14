      MODULE m_julia
      CONTAINS
      SUBROUTINE julia(
     >                 nop,nop2,mrot,tau,amat,bmat,omtil,jspins,latnam,
     >                 film,tria,invs,invs2,l_ss,l_soc,ndir,name,
     X                 nkpt,nmop,l_q)
!----------------------------------------------------------------------+
! Generate a k-point file with approx. nkpt k-pts or a Monkhorst-Pack  |
! set with nmod(i) divisions in i=x,y,z direction. Interface to kptmop |
! and kpttet routines of the MD-programm.                              |
!                                                          G.B. 07/01  |
!----------------------------------------------------------------------+

      USE m_constants, ONLY : pimach
      USE m_bravais
      USE m_divi
      USE m_brzone
      USE m_kptmop
      USE m_kpttet
      USE m_bandstr1
      IMPLICIT NONE

      INTEGER, PARAMETER :: nop48  = 48
      INTEGER, PARAMETER :: mface  = 51
      INTEGER, PARAMETER :: mdir   = 10
      INTEGER, PARAMETER :: nbsz   =  3
      INTEGER, PARAMETER :: ibfile =  6
      INTEGER, PARAMETER :: nv48   = (2*nbsz+1)**3+48

      INTEGER, INTENT (IN) :: nop,nop2,jspins,ndir
      REAL,    INTENT (IN) :: omtil
      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),tau(3,nop)
      LOGICAL, INTENT (IN) :: film,invs,invs2,l_ss,l_soc,l_q
      LOGICAL, INTENT (INOUT) :: tria
      CHARACTER(len=3), INTENT (IN) :: latnam
      CHARACTER(len=8), INTENT (IN) :: name(10)
      INTEGER, INTENT (INOUT) :: nkpt    ! number of k-points generated in set
      INTEGER, INTENT (INOUT) :: nmop(3) ! integer number triple: nmop(i), i=1,3;
                                         ! number of k-points in direction of rltv(ix,i)
      INTEGER ndiv3              ! max. number of tetrahedrons (< 6*(nkpt+1)
      INTEGER ntet               ! actual number of tetrahedrons 
      REAL, ALLOCATABLE :: vkxyz(:,:) ! vector of kpoint generated; in cartesian representation
      REAL, ALLOCATABLE ::  wghtkp(:) ! weight associated with k-points for BZ integration
      INTEGER, ALLOCATABLE :: ntetra(:,:) ! corners of the tetrahedrons
      REAL, ALLOCATABLE ::  voltet(:)     ! voulmes of the tetrahedrons
      REAL, ALLOCATABLE :: vktet(:,:)     !

      REAL    divis(4)           ! Used to find more accurate representation of k-points
                                 ! vklmn(i,kpt)/divis(i) and weights as wght(kpt)/divis(4)
      INTEGER nkstar             ! number of stars for k-points generated in full stars
      INTEGER nsym               ! number of symmetry elements
      REAL    bltv(3,3)          ! cartesian Bravais lattice basis (a.u.)
      REAL    rltv(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
      REAL    ccr(3,3,nop48)     ! rotation matrices in cartesian repr.
      REAL    rlsymr(3,3,nop48)  ! rotation matrices in reciprocal lattice basis representation
      REAL    talfa(3,nop48)     ! translation vector associated with (non-symmorphic)
                                 ! symmetry elements in Bravais lattice representation
      INTEGER ncorn,nedge,nface  ! number of corners, faces and edges of the IBZ
      REAL    fnorm(3,mface)     ! normal vector of the planes bordering the IBZ
      REAL    fdist(mface)       ! distance vector of the planes bordering t IBZ
      REAL    cpoint(3,mface)    ! cartesian coordinates of corner points of IBZ
      REAL    xvec(3)            ! arbitrary vector lying in the IBZ

      INTEGER idsyst   ! crystal system identification in MDDFT programs
      INTEGER idtype   ! lattice type identification in MDDFT programs

      INTEGER idimens  ! number of dimensions for k-point set (2 or 3)
      INTEGER nreg     ! 1 kpoints in full BZ; 0 kpoints in irrBZ
      INTEGER nfulst   ! 1 kpoints ordered in full stars
                       !    (meaningful only for nreg =1; full BZ)
      INTEGER nbound   ! 0 no primary points on BZ boundary;
                       ! 1 with boundary points (not for BZ integration!!!)
      INTEGER ikzero   ! 0 no shift of k-points;
                       ! 1 shift of k-points for better use of sym in irrBZ
      REAL    kzero(3)           ! shifting vector to bring one k-point to or 
                                 ! away from (0,0,0) (for even/odd nmop)

      INTEGER i,j,k,l,idiv,mkpt
      INTEGER iofile,iokpt,kpri,ktest,kmidtet
      INTEGER idivis(3)
      LOGICAL random,trias
      REAL help(3),tpi,binv(3,3),rlsymr1(3,3),ccr1(3,3)

      tpi = 2.0*pimach()
      random = .false.  ! do not use random tetra-points

c------------------------------------------------------------
c
c        idsyst         idtype 
c
c   1  cubic          primitive
c   2  tetragonal     body centered
c   3  orthorhombic   face centered
c   4  hexagonal      A-face centered
c   5  trigonal       B-face centered
c   6  monoclinic     C-face centered
c   7  triclinic 
c
c --->   for 2 dimensions only the following Bravais lattices exist:
c
c    TYPE                    EQUIVALENT 3-DIM        idsyst/idtype
c   square               = p-tetragonal ( 1+2 axis )      2/1
c   rectangular          = p-orthorhomb ( 1+2 axis )      3/1
c   centered rectangular = c-face-orthorhomb( 1+2 axis)   3/6
c   hexagonal            = p-hexagonal  ( 1+2 axis )      4/1
c   oblique              = p-monoclinic ( 1+2 axis )      6/1
c
c------------------------------------------------------------
      IF(l_q) THEN
       trias=tria
       tria=.false.
      ENDIF
       
      IF (latnam.EQ.'squ') THEN
        idsyst = 2
        idtype = 1
        IF (.not.film) THEN
          IF (abs(amat(1,1)-amat(3,3)) < 0.0000001) THEN
            idsyst = 1
            idtype = 1
          ENDIF
        ENDIF
      END IF
      IF (latnam.EQ.'p-r') THEN
        idsyst = 3
        idtype = 1
      END IF
      IF ((latnam.EQ.'c-b').OR.(latnam.EQ.'c-r')) THEN
        idsyst = 3
        idtype = 6
      END IF
      IF ((latnam.EQ.'hex').OR.(latnam.EQ.'hx3')) THEN
        idsyst = 4
        idtype = 1
      END IF
      IF (latnam.EQ.'obl') THEN
        idsyst = 6
        idtype = 1
      END IF
      IF (latnam.EQ.'any') THEN
        CALL bravais(
     >               amat,
     <               idsyst,idtype)
      ENDIF
      nsym = nop
      IF (film) nsym = nop2        
c
c-------------------- Want to make a Bandstructure ? --------
c
      IF (ndir == -4) THEN
        CALL bandstr1(
     >                idsyst,idtype,bmat,nkpt,name,jspins,film)
        RETURN
      ENDIF
c
c-------------------- Some variables we do not use ----------
c
      iofile = 6
      iokpt  = 6
      kpri   = 0 ! 3
      ktest  = 0 ! 5
      kmidtet = 0
      nreg    = 0
      nfulst  = 0
      ikzero  = 0
      kzero(1) = 0.0 ; kzero(2) = 0.0 ; kzero(3) = 0.0 
      nbound  = 0
      IF (tria) THEN
        IF (film) nbound  = 1
!        IF ((idsyst==1).AND.(idtype==1)) nbound  = 1
!        IF ((idsyst==2).AND.(idtype==1)) nbound  = 1
!        IF ((idsyst==3).AND.(idtype==1)) nbound  = 1
!        IF ((idsyst==3).AND.(idtype==6)) nbound  = 1
!        IF ((idsyst==4).AND.(idtype==1)) nbound  = 1
        IF (nbound == 0) random = .true.
      ENDIF
      idimens = 3
      IF (film) idimens = 2
c
c--------------------- Lattice information ------------------

      DO j = 1,3
        DO k = 1,3
          bltv(j,k) = amat(k,j)
          binv(j,k) = bmat(k,j)/tpi
          rltv(j,k) = bmat(k,j)
          DO i = 1,nsym
            rlsymr(k,j,i) = real( mrot(j,k,i) )
          ENDDO
        ENDDO
      ENDDO

      DO i = 1,nsym
        DO j = 1,3
          talfa(j,i) = 0.0
          DO k = 1,3
            talfa(j,i) = bltv(j,k) * tau(k,i)
            help(k) = 0.0
            DO l = 1,3
              help(k) =  help(k) + rlsymr(l,k,i) * binv(j,l)
            ENDDO
          ENDDO
          DO k = 1,3
           ccr(j,k,i) = 0.0
           DO l = 1,3
              ccr(j,k,i) = ccr(j,k,i) + bltv(l,k) * help(l)
            ENDDO
          ENDDO
        ENDDO
c      write (*,'(3f12.6)') ((ccr(j,k,i),j=1,3),k=1,3)
c      write (*,*)
      ENDDO
      DO i = 1,nsym
        rlsymr1(:,:) = rlsymr(:,:,i)
           ccr1(:,:) =    ccr(:,:,i)
        DO j = 1,3
          DO k = 1,3
            rlsymr(k,j,i) = rlsymr1(j,k)
               ccr(k,j,i) =    ccr1(j,k)
          ENDDO
        ENDDO
      ENDDO

      IF ((.not.l_ss).AND.(.not.l_soc).AND.(2*nsym<nop48)) THEN

        IF ( (film.AND.(.not.invs2)).OR.
     +     ((.not.film).AND.(.not.invs)) ) THEN
           ccr(:,:,nsym+1:2*nsym )    = -ccr(:,:,1:nsym)
           rlsymr(:,:,nsym+1:2*nsym ) = -rlsymr(:,:,1:nsym)
           nsym = 2 * nsym
        ENDIF

      ENDIF

C
C This subroutine finds the corner-points, the edges, and the
C faces of the irreducible wedge of the brillouin zone (IBZ).
C
      CALL brzone(
     >            rltv,nsym,ccr,mface,nbsz,nv48,
     =            cpoint,
     <            xvec,ncorn,nedge,nface,fnorm,fdist)

      IF ( tria.AND.random ) THEN
c
c       Calculate the points for tetrahedron method
c      
        mkpt = nkpt
        ndiv3 = 6*(mkpt+1)
        ALLOCATE (vkxyz(3,mkpt),wghtkp(mkpt) )
        ALLOCATE ( voltet(ndiv3),vktet(3,mkpt),ntetra(4,ndiv3) )
        CALL kpttet(
     >              iofile,ibfile,iokpt,
     >              kpri,ktest,kmidtet,mkpt,ndiv3,
     >              nreg,nfulst,rltv,omtil,
     >              nsym,ccr,mdir,mface,
     >              ncorn,nface,fdist,fnorm,cpoint,
     <              voltet,ntetra,ntet,vktet,
     =              nkpt,
     <              divis,vkxyz,wghtkp)
      ELSE
c
c       If just the total number of k-points is given, determine 
c       the divisions in each direction (nmop):
c
!        IF (tria) THEN
!            nkpt = nkpt/4
!            nmop(:) = nmop(:) / 2
!        ENDIF
        IF (nmop(1)+nmop(2)+nmop(3).EQ.0) THEN
          CALL divi(
     >              nkpt,bmat,film,nop,nop2,
     <              nmop)
        ENDIF
c
c       Now calculate Monkhorst-Pack k-points:
c
        IF (nmop(2).EQ.0) nmop(2) = nmop(1)
        IF ((.not.film).AND.(nmop(3).EQ.0)) nmop(3) = nmop(2)
        IF (nbound.EQ.1) THEN
           mkpt = (2*nmop(1)+1)*(2*nmop(2)+1)
           IF (.not.film) mkpt = mkpt*(2*nmop(3)+1)
        ELSE
           mkpt = nmop(1)*nmop(2)
           IF (.not.film) mkpt = mkpt*nmop(3)
        ENDIF
        ALLOCATE (vkxyz(3,mkpt),wghtkp(mkpt) )

        CALL kptmop(
     >              iofile,iokpt,kpri,ktest,
     >              idsyst,idtype,nmop,ikzero,kzero,
     >              rltv,bltv,nreg,nfulst,nbound,idimens,
     >              xvec,fnorm,fdist,ncorn,nface,nedge,cpoint,
     >              nsym,ccr,rlsymr,talfa,mkpt,mface,mdir,
     <              nkpt,divis,vkxyz,nkstar,wghtkp)

      ENDIF
c
      idivis(1) = int(divis(1)) 
      idivis(2) = int(divis(2)) 
      idivis(3) = int(divis(3)) 
      idiv = lcm(3,idivis)
c      WRITE (*,'(2i5)') nkpt,idiv
      IF (idiv.GE.200) idiv = 1
      DO j=1,nkpt
c        WRITE (*,'(4f10.5)') (vkxyz(i,j),i=1,3),wghtkp(j)
        wghtkp(j) = wghtkp(j) * divis(4)
        DO k = 1,3
          help(k) = 0.0
          DO l = 1,3
             help(k) = help(k) + amat(l,k) * vkxyz(l,j)
          ENDDO
        ENDDO
        DO i=1,3
          vkxyz(i,j) = help(i) * idiv / tpi
        ENDDO
      ENDDO
c
c if (l_q) write qpts file:
c
      IF(l_q)THEN
        IF(film) STOP'For the case of film q-points generator not
     &  implemented!'
        OPEN(113,file='qpts',form='formatted',status='new')
        WRITE(113,'(i5)') nkpt+1
        WRITE(113,8050) 0.,0.,0.
        DO j = 1, nkpt
           WRITE (113,FMT=8050) (vkxyz(i,j)/real(idiv),i=1,3)
        ENDDO
        CLOSE(113)
        tria=trias
        RETURN
      ENDIF
 8050 FORMAT (2(f14.10,1x),f14.10)

c
c write k-points file
c
      OPEN (41,file='kpts',form='formatted',status='new')
      IF (film) THEN
        WRITE (41,FMT=8110) nkpt,real(idiv),.false.
        DO j=nkpt,1,-1
           WRITE (41,FMT=8040) (vkxyz(i,j),i=1,2),wghtkp(j)
        ENDDO
      ELSE
        WRITE (41,FMT=8100) nkpt,real(idiv)
        DO j = 1, nkpt
           WRITE (41,FMT=8040) (vkxyz(i,j),i=1,3),wghtkp(j)
        ENDDO
        IF (tria.AND.random) THEN
          WRITE (41,'(i5)') ntet
          WRITE (41,'(4(4i6,4x))') ((ntetra(i,j),i=1,4),j=1,ntet)
          WRITE (41,'(4f20.13)') (ABS(voltet(j)),j=1,ntet)
        ENDIF
      END IF
 8100 FORMAT (i5,f20.10)
 8110 FORMAT (i5,f20.10,3x,l1)
 8040 FORMAT (4f10.5)
      CLOSE (41)

      DEALLOCATE ( vkxyz,wghtkp )
      IF (tria.AND..not.film)  DEALLOCATE ( voltet,vktet,ntetra )
      RETURN

      CONTAINS

      INTEGER FUNCTION lcm( n, ints )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute least common multiple (lcm) of n positive integers.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!===> Arguments
      INTEGER :: n
      INTEGER :: ints(n)

!===> Variables
      INTEGER :: i,j,m
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF ( any( ints(1:n)<= 0 ) ) THEN
        m = 0
      ELSE
        m = maxval( ints(1:n) )
        DO i = 1, n
          DO j = 1, ints(i)/2
            IF ( mod( m*j,ints(i) ) == 0 ) EXIT
          END DO
          m = m*j
        ENDDO
      ENDIF

      lcm = m

      RETURN
      END FUNCTION lcm

      END SUBROUTINE julia
      END MODULE m_julia
