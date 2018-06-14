      MODULE m_set
      CONTAINS
!---------------------------------------------- 
      SUBROUTINE cset(idim,const,array,istride)

      IMPLICIT NONE

      INTEGER idim,istride
      COMPLEX const,array(idim)

      INTEGER i

      DO i=1,idim,istride
         array(i)=const
      ENDDO

      END SUBROUTINE cset
!---------------------------------------------- 
      SUBROUTINE dset(idim,const,array,istride)

      IMPLICIT NONE

      INTEGER idim,istride
      REAL    const,array(idim)

      INTEGER i

      DO i=1,idim,istride
         array(i)=const
      ENDDO

      END SUBROUTINE dset
!---------------------------------------------- 
      END MODULE m_set
