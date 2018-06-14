      MODULE m_xsf_io
!-----------------------------------------------
! DESC:subroutines to write xsf-files for xcrysden
!                 Daniel Wortmann, (06-01-26)
!-----------------------------------------------
      CONTAINS
      !<-- S:S: xsf_WRITE_atoms(fileno,film,amat,neq(:ntype),zatom(:ntype),pos)
      SUBROUTINE xsf_WRITE_atoms(fileno,film,od,amat,neq,zatom,pos)
!-----------------------------------------------
!     Writes the crystal dimensions&atomic positions
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: fileno
      LOGICAL,INTENT(IN)     :: film,od
      REAL   ,INTENT(IN)     :: amat(3,3)
      INTEGER,INTENT(IN)     :: neq(:)
      REAL   ,INTENT(IN)     :: zatom(:)
      REAL   ,INTENT(IN)     :: pos(:,:)
      !>
      !<-- Locals
      INTEGER             :: n,nn,na
      !>
      IF (film) THEN
         IF (od) THEN
            WRITE(fileno,*) "POLYMERE"
         ELSE
            WRITE(fileno,*) "SLAB"
         ENDIF
      ELSE
         WRITE(fileno,*) "CRYSTAL"
      ENDIF

      WRITE(fileno,*) "PRIMVEC"
      WRITE(fileno,'(3(f0.7,1x))') amat(:,1)
      WRITE(fileno,'(3(f0.7,1x))') amat(:,2)
      WRITE(fileno,'(3(f0.7,1x))') amat(:,3)
      
      WRITE(fileno,*) "PRIMCOORD"
      WRITE(fileno,*) SUM(neq)," 1"
      na = 1
      DO n = 1,SIZE(neq)
         DO nn = 1,neq(n)
            WRITE(fileno,'(i4,2x,3(f0.7,1x))') NINT(zatom(n)),pos(:,na)
            na=na+1
         ENDDO
      ENDDO
      WRITE(fileno,*)
      END SUBROUTINE
      !> 
      !<-- S: xsf_write_header(fileno,twodim,desc,vec1,vec2,vec3,zero,grid)
      SUBROUTINE xsf_WRITE_header(fileno,twodim,desc,vec1,vec2,vec3,zero
     $     ,grid)
!-----------------------------------------------
!  writes the beginning of a gid-datablock
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: fileno,grid(:)
      LOGICAL,INTENT(IN)     :: twodim
      REAL   ,INTENT(IN)     :: vec1(:),vec2(:),vec3(:),zero(:)
      CHARACTER(LEN =*),INTENT(IN) :: desc 
      !>

      if (twodim) THEN
         WRITE(fileno,*) "BEGIN_BLOCK_DATAGRID_2D"
         WRITE(fileno,*) desc
         WRITE(fileno,*) "BEGIN_DATAGRID_2D_A"
         WRITE(fileno,'(3i7)') grid(1:2)
         WRITE(fileno,'(3(f12.7,1x))') zero
         WRITE(fileno,'(3(f12.7,1x))') vec1
         WRITE(fileno,'(3(f12.7,1x))') vec2 
      ELSE
         WRITE(fileno,*) "BEGIN_BLOCK_DATAGRID_3D"
         WRITE(fileno,*) desc
         WRITE(fileno,*) "BEGIN_DATAGRID_3D_A"
         WRITE(fileno,'(3i7)') grid(1:3)
         WRITE(fileno,'(3(f12.7,1x))') zero
         WRITE(fileno,'(3(f12.7,1x))') vec1
         WRITE(fileno,'(3(f12.7,1x))') vec2
         WRITE(fileno,'(3(f12.7,1x))') vec3    
      ENDIF
      END SUBROUTINE
      !> 
      !<-- S: xsf_write_newblock(fileno,twodim,vec1,vec2,vec3,zero,grid)
      SUBROUTINE xsf_WRITE_newblock(fileno,twodim,vec1,vec2
     $     ,vec3,zero,grid)
!-----------------------------------------------
!  writes the beginning of a new gid-datablock for second spin
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: fileno,grid(:)
      LOGICAL,INTENT(IN)     :: twodim
      REAL   ,INTENT(IN)     :: vec1(:),vec2(:),vec3(:),zero(:)
      !>

      if (twodim) THEN
         WRITE(fileno,*) "END_DATAGRID_2D"
         WRITE(fileno,*) "BEGIN_DATAGRID_2D_B"
         WRITE(fileno,'(3i7)') grid(1:2)
         WRITE(fileno,'(3(f12.7,1x))') zero
         WRITE(fileno,'(3(f12.7,1x))') vec1
         WRITE(fileno,'(3(f12.7,1x))') vec2
      ELSE
         WRITE(fileno,*) "END_DATAGRID_3D"
         WRITE(fileno,*) "BEGIN_DATAGRID_3D_B"
         WRITE(fileno,'(3i7)') grid(1:3)
         WRITE(fileno,'(3(f12.7,1x))') zero
         WRITE(fileno,'(3(f12.7,1x))') vec1
         WRITE(fileno,'(3(f12.7,1x))') vec2
         WRITE(fileno,'(3(f12.7,1x))') vec3
      ENDIF
      END SUBROUTINE
      !> 
      !<-- S: xsf_write_endblock(fileno,twodim)
      SUBROUTINE xsf_write_endblock(fileno,twodim)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: fileno
      LOGICAL,INTENT(IN)     :: twodim
      !>
 
      IF (twodim) THEN
         WRITE(fileno,*) "END_DATAGRID_2D"                      
         WRITE(fileno,*) "END_BLOCK_DATAGRID_2D" 
      ELSE
         WRITE(fileno,*) "END_DATAGRID_3D"                      
         WRITE(fileno,*) "END_BLOCK_DATAGRID_3D" 
      ENDIF
      END SUBROUTINE
      !> 
 
      SUBROUTINE xsf_WRITE_force(fileno,film,od,amat,neq,zatom,pos,
     >                                                       force)
!-----------------------------------------------
!     Writes the crystal dimensions&force positions
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: fileno
      LOGICAL,INTENT(IN)     :: film,od
      REAL   ,INTENT(IN)     :: amat(3,3)
      INTEGER,INTENT(IN)     :: neq(:)
      REAL   ,INTENT(IN)     :: zatom(:)
      REAL   ,INTENT(IN)     :: pos(:,:)
      INTEGER,INTENT(IN)     :: force ! number of atoms + force vectors
      !>
      !<-- Locals
      INTEGER             :: n,nn,na
      !>
      IF (film) THEN
         IF (od) THEN
            WRITE(fileno,*) "POLYMERE"
         ELSE
            WRITE(fileno,*) "SLAB"
         ENDIF
      ELSE
         WRITE(fileno,*) "CRYSTAL"
      ENDIF

      WRITE(fileno,*) "PRIMVEC"
      WRITE(fileno,'(3(f0.7,1x))') amat(:,1)
      WRITE(fileno,'(3(f0.7,1x))') amat(:,2)
      WRITE(fileno,'(3(f0.7,1x))') amat(:,3)
      
      WRITE(fileno,*) "PRIMCOORD"
      WRITE(fileno,*) force," 1"
      na = 1
      DO n = 1,SIZE(neq)
         DO nn = 1,neq(n)
            WRITE(fileno,'(i4,2x,3(f0.7,1x))') NINT(zatom(n)),pos(:,na)
            na=na+1
         ENDDO
      ENDDO
      WRITE(fileno,*)
      END SUBROUTINE
      !> 
!-----------------------------------------------
      END MODULE

 
