      PROGRAM inpgen
!----------------------------------------------------------------------------+
!   Set up a FLEUR inp-file from basic input data; for use and docu please   !
!   refer to inpgen.html (or see http://www.flapw.de/docs/inpgen.html)       !
!                                                                            |
!   The program is based on the input conventions of the FLAIR-code, so that !
!   some compatibility is ensured. The symmetry generator was written by     ! 
!   M.Weinert and implemented in the FLAIR-code by G.Schneider.              !
!                                                                    gb`02   |
!----------------------------------------------------------------------------+
  
      USE m_structinput
      USE m_crystal
      USE m_socorssdw
      USE m_rwsymfile
      USE m_setinp
      USE m_writestruct
      USE m_xsf_io, ONLY : xsf_write_atoms

      IMPLICIT NONE
    
      INTEGER natmax,nop48,nline,natin,ngen,i,j
      INTEGER nop,nop2,ntype,nat,nops,no3,no2,na
      INTEGER infh,errfh,bfh,warnfh,symfh,dbgfh,outfh,dispfh
      LOGICAL film,cal_symm,checkinp,symor,invs,zrfs,invs2
      LOGICAL cartesian,oldfleur,l_soc,l_ss,inistop
      REAL    aa,dvac,theta,phi,omtil
 
      REAL a1(3),a2(3),a3(3),scale(3),amat(3,3),bmat(3,3),qss(3)
      INTEGER, ALLOCATABLE :: mmrot(:,:,:)
      REAL,    ALLOCATABLE :: ttr(:,:),atompos(:,:),atomid(:),taual(:,:)
      REAL,    ALLOCATABLE :: rmt(:)
      INTEGER, POINTER :: neq(:), ntyrep(:)              ! these variables are allocated with
      REAL,    POINTER :: zatom(:)                       ! dim 'ntype'
      INTEGER, POINTER :: natype(:),natrep(:),natmap(:)  ! or  'nat'
      REAL,    POINTER :: pos(:,:)                       ! or  '3,nat'
      INTEGER, POINTER :: mrot(:,:,:)                    ! or  '3,3,nop'
      REAL,    POINTER :: tau(:,:)                       ! or  '3,nop' here, or in atom_sym

      CHARACTER(len=80):: title
      CHARACTER(len=7) :: symfn
      CHARACTER(len=4) :: dispfn

      nop48 = 48
      natmax = 999
      ngen = 0
      infh = 5
      errfh = 6 ; warnfh = 6 ; dbgfh = 6 ; outfh = 6
      bfh = 93
      symfh = 94
      symfn = 'sym    '
      dispfh = 97
      dispfn='disp'
      nline = 0

      ALLOCATE ( mmrot(3,3,nop48), ttr(3,nop48) )
      ALLOCATE ( atompos(3,natmax),atomid(natmax) )

!      OPEN (5,file='inp2',form='formatted',status='old')
      OPEN (6,file='out',form='formatted',status='unknown')

      CALL struct_input(
     >                  infh,errfh,bfh,warnfh,symfh,symfn,
     >                  natmax,nop48,
     X                  nline,
     <                  title,film,cal_symm,checkinp,symor,
     <                  cartesian,oldfleur,a1,a2,a3,dvac,aa,scale,
     <                  natin,atomid,atompos,ngen,mmrot,ttr,
     <                  l_soc,l_ss,theta,phi,qss,inistop)

!      CLOSE (5)

      IF (.not.film) dvac=a3(3)
      WRITE (6,*)
      WRITE (6,*) title
      WRITE (6,*) 'film=',film,'cal_symm=',cal_symm
      WRITE (6,*) 'checkinp=',checkinp,'symor=',symor
      WRITE (6,*)
      WRITE (6,'(a5,3f10.5)') 'a1 = ',a1(:)
      WRITE (6,'(a5,3f10.5)') 'a2 = ',a2(:)
      WRITE (6,'(a5,3f10.5)') 'a3 = ',a3(:)
      WRITE (6,*)
      WRITE (6,'(2(a5,f10.5))') 'dvac=',dvac,' aa =',aa
      WRITE (6,'(a8,3f10.5)') 'scale = ',scale(:)
      WRITE (6,*)
      WRITE (6,'(a6,i3,a6,999i5)') 'natin=',natin,' Z = ',
     +                             (nint(atomid(i)),i=1,abs(natin))
      WRITE (6,*) 'positions: '
      WRITE (6,'(3(3x,f10.5))') ((atompos(j,i),j=1,3),i=1,abs(natin))
      WRITE (6,*)
      WRITE (6,*) 'generators: ',ngen,'(excluding identity)'
      DO i = 2, ngen+1
         WRITE (6,*) i
         WRITE (6,'(3i5,f8.3)') (mmrot(1,j,i),j=1,3),ttr(1,i)
         WRITE (6,'(3i5,f8.3)') (mmrot(2,j,i),j=1,3),ttr(2,i)
         WRITE (6,'(3i5,f8.3)') (mmrot(3,j,i),j=1,3),ttr(3,i)
      ENDDO
      IF (l_soc) WRITE(6,'(a4,2f10.5)') 'soc:',theta,phi
      IF (l_ss)  WRITE(6,'(a4,3f10.5)') 'qss:',qss(:)
!
! --> generate symmetry from input (atomic positions, generators or whatever)
!     
      CALL crystal(
     >             dbgfh,errfh,outfh,dispfh,dispfn,
     >             cal_symm,cartesian,symor,film, ! gb: changed: oldfleur,
     >             natin,natmax,nop48,
     >             atomid,atompos,a1,a2,a3,aa,scale,
     <             invs,zrfs,invs2,nop,nop2,
     <             ngen,mmrot,ttr,ntype,nat,nops,
     <             neq,ntyrep,zatom,natype,natrep,natmap,
     <             mrot,tau,pos,amat,bmat,omtil)

      IF (l_ss.OR.l_soc)  THEN
         CALL soc_or_ssdw(
     >                    l_soc,l_ss,theta,phi,qss,amat,
     >                    mrot,tau,nop,nop2,nat,atomid,atompos,
     <                    mmrot,ttr,no3,no2,ntype,neq,natmap,
     <                    ntyrep,natype,natrep,zatom,pos)
         nop = no3 ; nop2 = no2
         mrot(:,:,1:nop) = mmrot(:,:,1:nop)
         tau(:,1:nop) = ttr(:,1:nop)
      ENDIF
      DEALLOCATE ( mmrot, ttr, atomid, atompos )

      ALLOCATE ( taual(3,nat) ) 
      WRITE (6,*)
      WRITE (6,'(a6,i3,a6,i3)') 'ntype=',ntype,' nat= ',nat
      na = 0
      DO i = 1, ntype
        WRITE (6,'(a3,i3,a2,i3,a6,i3)') ' Z(',i,')=',nint(zatom(i)),
     +                                             ' neq= ',neq(i)
        DO j = 1, neq(i)
           WRITE (6,'(3f10.6,10x,i7)')
     &           pos(:,natmap(na+j)),natmap(na+j)
           taual(:,na+j) = pos(:,natmap(na+j))
        ENDDO
        na = na + neq(i)
      ENDDO
      DO i=1,nat
        pos(:,i) = matmul( amat , taual(:,i) )
      ENDDO
      DEALLOCATE ( ntyrep, natype, natrep )
!
! --> write a file 'sym.out' with accepted symmetry operations
!
      nops = nop
      symfn = 'sym.out'
      IF (.not.film) nop2=nop
      CALL rw_symfile(
     >                'W',symfh,symfn,nops,bmat,
     X                 mrot,tau,nop,nop2,symor)

!
! --> set defaults for FLEUR inp-file
!
      ALLOCATE ( rmt(ntype) )
      CALL set_inp(
     >             ntype,nop,nop2,nat,omtil,title,
     >             film,invs,invs2,zrfs,amat,dvac,l_soc,
     >             neq,zatom,taual,l_ss,a1,a2,a3,theta,phi,
     <             rmt)
!
! --> Structure in povray or xsf-format
!
      IF (.false.) THEN
         CALL write_struct(
     >                  ntype,nat,neq,
     >                  rmt,pos,natmap,amat)
      ELSE 
         OPEN (55,file="struct.xsf")
         CALL xsf_WRITE_atoms(
     >                        55,film,.false.,amat,neq(:ntype),
     >                        zatom(:ntype),pos)
         CLOSE (55)
      ENDIF

      DEALLOCATE ( taual, mrot, tau, neq, zatom, rmt, natmap, pos )

      IF (inistop) STOP 'symmetry done' 

      END 
