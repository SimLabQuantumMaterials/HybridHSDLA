      MODULE m_constants
      IMPLICIT NONE
      CONTAINS
!------------------------------------------------------------------------
      REAL PURE FUNCTION pimach()
!
!     This subprogram supplies the value of the constant PI correct to
!     machine precision where
!
!     PI=3.1415926535897932384626433832795028841971693993751058209749446
!
      pimach = 3.1415926535897932
      END FUNCTION pimach
!------------------------------------------------------------------------
      REAL ELEMENTAL FUNCTION c_light(fac)
!
!     This subprogram supplies the value of c according to
!     NIST standard 13.1.99 
!     Hartree and Rydbergs changed by fac = 1.0 or 2.0
!
      REAL, INTENT (IN) :: fac
      c_light = 137.0359895e0 * fac 

      RETURN
      END FUNCTION c_light
      END MODULE m_constants

