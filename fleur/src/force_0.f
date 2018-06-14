      MODULE m_force0
      CONTAINS
      SUBROUTINE force_0(
     >                   jspd,ntypd,
     X                   force,
     <                   force_old)
c ************************************************************
c Presetting force components
c Also presets lm quantities for Pulay terms
c al la Yu et al equa A12,A17,A20
c ************************************************************
c
      IMPLICIT NONE 
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: jspd,ntypd
C     ..
C     .. Arrays Arguments ..
      REAL,   INTENT (INOUT):: force(3,ntypd,jspd)
      REAL,   INTENT (OUT)  :: force_old(3,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER i,jsp,n
      REAL zero
C     ..
C     .. Data statements ..
      DATA zero/0.0/
C     ..
c     preset force components to zero and save force
c
      DO n = 1,ntypd
         DO i = 1,3
           force_old(i,n) = zero
           DO jsp = 1,jspd
               force_old(i,n) = force_old(i,n) + force(i,n,jsp)
               force(i,n,jsp) = zero
            END DO
         END DO
      END DO
c
      END SUBROUTINE force_0
      END MODULE m_force0
