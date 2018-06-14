      MODULE m_setmn
c.....------------------------------------------------------------------
c     set minimum value to af2 if it is smaller than some small value.
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE setmn(
     >                 ifftd2,kmx,jspins,
     X                 af2)

      IMPLICIT NONE
      INTEGER, INTENT (IN)    :: ifftd2,jspins,kmx
      REAL,    INTENT (INOUT) :: af2(0:ifftd2-1,jspins)

      INTEGER i,js
      REAL, PARAMETER :: sml = 1.e-13
      INTRINSIC max

      DO js = 1,jspins
          DO i = 0,kmx
              af2(i,js) = max(af2(i,js),sml)
          ENDDO
      ENDDO

      END SUBROUTINE setmn
      END MODULE m_setmn
