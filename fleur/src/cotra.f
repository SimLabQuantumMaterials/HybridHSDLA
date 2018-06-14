      MODULE m_cotra

!******************************************************************
! transform, in real (reciprocal) space, internal coordinates r (f)
! into cartesian coordinates x (g) or back:
!
! cotra0   real-space  internal  -> cartesian    r -> x
! cotra1   real-space  cartesian -> internal     x -> r
! cotra3   reciprocal  internal  -> cartesian    f -> g
!******************************************************************

      IMPLICIT NONE
 
      REAL,    PRIVATE :: s
      INTEGER, PRIVATE :: i,j

      CONTAINS

!     **********************************************************
      SUBROUTINE cotra0(r,x,amat)
!     **********************************************************
!     ..
!     .. Array Arguments ..
      REAL, INTENT (IN) :: r(3)
      REAL, INTENT (IN) :: amat(3,3)
      REAL, INTENT (OUT):: x(3)
!
      DO i = 1,3
         s = 0.0
         DO j = 1,3
            s = s + amat(i,j)*r(j)
         ENDDO
         x(i) = s
      ENDDO

      RETURN
      END SUBROUTINE cotra0

!     ***********************************************************
      SUBROUTINE cotra1(x,r,bmat)
!     ***********************************************************
      USE m_constants, ONLY : pimach
!     ..
!     .. Array Arguments ..
      REAL, INTENT (IN) :: x(3)
      REAL, INTENT (IN) :: bmat(3,3)
      REAL, INTENT (OUT):: r(3)
!     ..
!     .. Local Scalars ..
      REAL tpi
!
      tpi = 2 * pimach()
!
      DO i = 1,3
         s = 0.0
         DO j = 1,3
            s = s + bmat(i,j)*x(j)
         ENDDO
         r(i) = s/tpi
      ENDDO

      RETURN
      END SUBROUTINE cotra1

!*********************************************************************
      SUBROUTINE cotra3(f,g,bmat)
!*********************************************************************
!     ..
!     .. Array Arguments ..
      REAL,    INTENT (IN) :: f(3)
      REAL,    INTENT (IN) :: bmat(3,3)
      REAL,    INTENT (OUT):: g(3)
!
      DO j = 1,3
         g(j) = 0.0
         DO i = 1,3
            g(j) = g(j) + f(i)*bmat(i,j)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE cotra3

      END MODULE m_cotra
