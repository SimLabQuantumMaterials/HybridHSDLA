      MODULE m_spgrot
!     **********************************************************
!     perform space group operations of film
!     **********************************************************
      CONTAINS
      SUBROUTINE spgrot(
     >                  nop,symor,tpi,mrot,tau,invtab,
     >                  k,
     <                  kr,phas)

      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: nop
      REAL,    INTENT (IN)  :: tpi
      LOGICAL, INTENT (IN)  :: symor
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: k(3),mrot(3,3,nop),invtab(nop)
      REAL,    INTENT (IN)  :: tau(3,nop)
      INTEGER, INTENT (OUT) :: kr(3,nop)
      REAL,    INTENT (OUT) :: phas(nop)
!     ..
!     .. Local Scalars ..
      INTEGER n,ni
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC cos
!     ..
      DO n = 1,nop
         kr(1,n) = k(1)*mrot(1,1,n) + k(2)*mrot(2,1,n) +
     +             k(3)*mrot(3,1,n)
         kr(2,n) = k(1)*mrot(1,2,n) + k(2)*mrot(2,2,n) +
     +             k(3)*mrot(3,2,n)
         kr(3,n) = k(1)*mrot(1,3,n) + k(2)*mrot(2,3,n) +
     +             k(3)*mrot(3,3,n)
      ENDDO
      IF (symor) THEN
         DO n = 1,nop
            phas(n) = 1.
         ENDDO
      ELSE
         DO n = 1,nop
            ni = invtab(n)
            phas(n) = cos(tpi* ( kr(1,n)*tau(1,ni) +
     +                           kr(2,n)*tau(2,ni) +
     +                           kr(3,n)*tau(3,ni)))
! note that, in general phas(n) could be complex!
         ENDDO
      END IF
      RETURN
      END SUBROUTINE spgrot
      END MODULE m_spgrot
