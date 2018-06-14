      SUBROUTINE cset(idim,const,array,istride)

      IMPLICIT NONE

      INTEGER idim,istride
      COMPLEX const,array(idim)

      INTEGER i

      DO i=1,idim,istride
         array(i)=const
      ENDDO

      RETURN
      END
